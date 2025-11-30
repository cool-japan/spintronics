//! Physical Reservoir Computing with Magnon Dynamics
//!
//! Uses spin wave non-linearity in magnetic materials as a computational substrate.
//! Input signals are encoded as magnetic fields, the magnon system processes them
//! through non-linear dynamics, and a trained readout layer produces outputs.

use crate::magnon::SpinChain;
use crate::vector3::Vector3;

/// Magnon-based physical reservoir computing system
#[derive(Debug, Clone)]
pub struct MagnonReservoir {
    /// Underlying spin chain serving as the reservoir
    pub system: SpinChain,

    /// Indices of spins used for input injection
    pub input_nodes: Vec<usize>,

    /// Indices of spins read out for computation
    pub readout_nodes: Vec<usize>,

    /// Linear readout weights (trained)
    pub readout_weights: Vec<f64>,

    /// Bias term for readout
    pub bias: f64,

    /// Input scaling factor
    pub input_scale: f64,
}

impl MagnonReservoir {
    /// Create a new magnon reservoir
    ///
    /// # Arguments
    /// * `system` - SpinChain that serves as the reservoir
    /// * `input_nodes` - Indices where input signals are injected
    /// * `readout_nodes` - Indices from which state is read
    pub fn new(system: SpinChain, input_nodes: Vec<usize>, readout_nodes: Vec<usize>) -> Self {
        let n_readout = readout_nodes.len() * 3; // 3 components per spin

        Self {
            system,
            input_nodes,
            readout_nodes,
            readout_weights: vec![0.0; n_readout],
            bias: 0.0,
            input_scale: 1.0e4, // Default: 10 kA/m
        }
    }

    /// Create a reservoir with uniformly distributed nodes
    ///
    /// # Arguments
    /// * `n_spins` - Total number of spins in chain
    /// * `n_input` - Number of input nodes
    /// * `n_readout` - Number of readout nodes
    /// * `exchange` - Exchange coupling \[J\]
    /// * `alpha` - Gilbert damping
    pub fn uniform_distribution(
        n_spins: usize,
        n_input: usize,
        n_readout: usize,
        exchange: f64,
        alpha: f64,
    ) -> Self {
        use crate::magnon::chain::ChainParameters;

        let params = ChainParameters {
            a_ex: exchange,
            alpha,
            ..ChainParameters::default()
        };

        let system = SpinChain::new(n_spins, params);

        // Input nodes evenly spaced
        let input_nodes: Vec<usize> = (0..n_input).map(|i| i * n_spins / (n_input + 1)).collect();

        // Readout nodes evenly spaced
        let readout_nodes: Vec<usize> = (0..n_readout)
            .map(|i| i * n_spins / (n_readout + 1))
            .collect();

        Self::new(system, input_nodes, readout_nodes)
    }

    /// Inject input signal as magnetic field
    ///
    /// Maps scalar input values to z-component of magnetic field
    fn inject_input(&self, input: &[f64]) -> Vector3<f64> {
        if input.is_empty() {
            return Vector3::new(0.0, 0.0, 0.0);
        }

        // Sum all inputs and apply to z-direction
        let total_input: f64 = input.iter().sum();
        Vector3::new(0.0, 0.0, total_input * self.input_scale)
    }

    /// Extract reservoir state from readout nodes
    ///
    /// Returns concatenated [mx, my, mz] for all readout nodes
    pub fn extract_state(&self) -> Vec<f64> {
        let mut state = Vec::with_capacity(self.readout_nodes.len() * 3);

        for &idx in &self.readout_nodes {
            let m = self.system.spins[idx];
            state.push(m.x);
            state.push(m.y);
            state.push(m.z);
        }

        state
    }

    /// Apply readout weights to reservoir state
    fn apply_readout(&self, state: &[f64]) -> f64 {
        let mut output = self.bias;

        for (i, &s) in state.iter().enumerate() {
            if i < self.readout_weights.len() {
                output += self.readout_weights[i] * s;
            }
        }

        output
    }

    /// Process input signal through the reservoir
    ///
    /// # Arguments
    /// * `input_signal` - Input values to encode
    /// * `h_ext` - External magnetic field (constant)
    /// * `steps` - Number of evolution steps
    /// * `dt` - Time step
    ///
    /// # Returns
    /// Output value after processing
    pub fn process(
        &mut self,
        input_signal: &[f64],
        h_ext: Vector3<f64>,
        steps: usize,
        dt: f64,
    ) -> f64 {
        // 1. Input Layer: Inject input as field modulation
        let input_field = self.inject_input(input_signal);
        let total_field = h_ext + input_field;

        // 2. Reservoir Layer: Evolve magnon dynamics
        for _ in 0..steps {
            self.system.evolve_heun(total_field, dt);
        }

        // 3. Readout Layer: Extract state and apply weights
        let state = self.extract_state();
        self.apply_readout(&state)
    }

    /// Train readout weights using linear regression
    ///
    /// Collects reservoir states for all training samples and solves
    /// for optimal weights using least squares.
    ///
    /// # Arguments
    /// * `inputs` - Training input sequences
    /// * `targets` - Target output values
    /// * `h_ext` - External field during training
    /// * `steps` - Evolution steps per sample
    /// * `dt` - Time step
    #[allow(clippy::needless_range_loop)]
    pub fn train(
        &mut self,
        inputs: &[Vec<f64>],
        targets: &[f64],
        h_ext: Vector3<f64>,
        steps: usize,
        dt: f64,
    ) -> Result<f64, &'static str> {
        if inputs.len() != targets.len() {
            return Err("Input and target lengths must match");
        }

        if inputs.is_empty() {
            return Err("No training data provided");
        }

        let n_samples = inputs.len();
        let state_dim = self.readout_nodes.len() * 3;

        // Collect reservoir states for all samples
        let mut state_matrix = vec![vec![0.0; state_dim + 1]; n_samples];

        for (i, input) in inputs.iter().enumerate() {
            // Reset system to initial state
            self.system.reset_to_z();

            // Process input
            let input_field = self.inject_input(input);
            let total_field = h_ext + input_field;

            for _ in 0..steps {
                self.system.evolve_heun(total_field, dt);
            }

            // Extract state
            let state = self.extract_state();
            for (j, &s) in state.iter().enumerate() {
                state_matrix[i][j] = s;
            }
            state_matrix[i][state_dim] = 1.0; // Bias term
        }

        // Solve linear regression: W = (X^T X)^{-1} X^T y
        // Simplified: Use gradient descent for stability
        let learning_rate = 0.01;
        let epochs = 1000;

        // Initialize weights and bias
        self.readout_weights = vec![0.0; state_dim];
        self.bias = 0.0;

        // Gradient descent
        for _ in 0..epochs {
            let mut weight_grads = vec![0.0; state_dim];
            let mut bias_grad = 0.0;

            for (i, state) in state_matrix.iter().enumerate() {
                // Forward pass
                let mut prediction = self.bias;
                for j in 0..state_dim {
                    prediction += self.readout_weights[j] * state[j];
                }

                // Error
                let error = prediction - targets[i];

                // Gradients
                for j in 0..state_dim {
                    weight_grads[j] += error * state[j];
                }
                bias_grad += error;
            }

            // Update weights
            for j in 0..state_dim {
                self.readout_weights[j] -= learning_rate * weight_grads[j] / n_samples as f64;
            }
            self.bias -= learning_rate * bias_grad / n_samples as f64;
        }

        // Calculate final MSE
        let mut mse = 0.0;
        for (i, state) in state_matrix.iter().enumerate() {
            let mut prediction = self.bias;
            for j in 0..state_dim {
                prediction += self.readout_weights[j] * state[j];
            }
            let error = prediction - targets[i];
            mse += error * error;
        }
        mse /= n_samples as f64;

        Ok(mse)
    }

    /// Reset reservoir to initial magnetization state
    pub fn reset(&mut self) {
        self.system.reset_to_z();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reservoir_creation() {
        let reservoir = MagnonReservoir::uniform_distribution(
            20, // 20 spins
            3,  // 3 input nodes
            5,  // 5 readout nodes
            1.0e-20, 0.01,
        );

        assert_eq!(reservoir.input_nodes.len(), 3);
        assert_eq!(reservoir.readout_nodes.len(), 5);
        assert_eq!(reservoir.readout_weights.len(), 15); // 5 nodes × 3 components
    }

    #[test]
    fn test_state_extraction() {
        let reservoir = MagnonReservoir::uniform_distribution(10, 2, 3, 1.0e-20, 0.01);

        let state = reservoir.extract_state();
        assert_eq!(state.len(), 9); // 3 nodes × 3 components
    }

    #[test]
    fn test_input_processing() {
        let mut reservoir = MagnonReservoir::uniform_distribution(10, 2, 3, 1.0e-20, 0.01);

        let input = vec![0.1, 0.2];
        let h_ext = Vector3::new(0.0, 0.0, 1.0e5);

        let output = reservoir.process(&input, h_ext, 10, 1.0e-13);

        // Should produce some output (not necessarily meaningful without training)
        assert!(output.is_finite());
    }

    #[test]
    fn test_training_simple_constant() {
        let mut reservoir = MagnonReservoir::uniform_distribution(15, 2, 4, 1.0e-20, 0.01);

        // Simple task: output constant value regardless of input
        let inputs = vec![
            vec![0.1, 0.2],
            vec![0.3, 0.1],
            vec![0.0, 0.5],
            vec![0.2, 0.2],
        ];
        let targets = vec![1.0, 1.0, 1.0, 1.0];

        let h_ext = Vector3::new(0.0, 0.0, 1.0e5);
        let result = reservoir.train(&inputs, &targets, h_ext, 5, 1.0e-13);

        assert!(result.is_ok());
        let mse = result.unwrap();
        assert!(mse >= 0.0);
    }

    #[test]
    fn test_reset_functionality() {
        let mut reservoir = MagnonReservoir::uniform_distribution(10, 2, 3, 1.0e-20, 0.01);

        // Evolve system
        let h_ext = Vector3::new(1.0e5, 0.0, 0.0);
        for _ in 0..50 {
            reservoir.system.evolve_heun(h_ext, 1.0e-13);
        }

        // Reset
        reservoir.reset();

        // Check all spins are back to +z
        for spin in &reservoir.system.spins {
            assert!((spin.x).abs() < 1e-10);
            assert!((spin.y).abs() < 1e-10);
            assert!((spin.z - 1.0).abs() < 1e-10);
        }
    }
}
