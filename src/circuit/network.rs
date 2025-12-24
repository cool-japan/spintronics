//! Spin circuit networks
//!
//! Multi-terminal spin circuits with Kirchhoff's laws for spin currents.

use std::collections::HashMap;

use super::resistor::SpinResistor;

/// Circuit terminal (node)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Terminal(pub usize);

/// Connection between terminals
#[derive(Debug, Clone)]
pub struct Connection {
    /// Terminal A
    pub from: Terminal,
    /// Terminal B
    pub to: Terminal,
    /// Resistor element
    pub resistor: SpinResistor,
}

/// Spin circuit network
///
/// Implements Kirchhoff's current laws for both charge and spin currents:
/// - ΣI_charge = 0 (charge conservation)
/// - ΣI_spin = -μ_s/τ_s (spin relaxation)
#[derive(Debug, Clone)]
pub struct SpinCircuit {
    /// Connections between terminals
    pub connections: Vec<Connection>,

    /// Applied voltages at terminals \[V\]
    pub voltages: HashMap<Terminal, f64>,

    /// Spin accumulation at terminals \[J\]
    pub spin_accumulation: HashMap<Terminal, f64>,
}

impl SpinCircuit {
    /// Create a new empty circuit
    pub fn new() -> Self {
        Self {
            connections: Vec::new(),
            voltages: HashMap::new(),
            spin_accumulation: HashMap::new(),
        }
    }

    /// Add a connection between two terminals
    pub fn add_connection(&mut self, from: Terminal, to: Terminal, resistor: SpinResistor) {
        self.connections.push(Connection { from, to, resistor });
    }

    /// Set voltage at a terminal
    pub fn set_voltage(&mut self, terminal: Terminal, voltage: f64) {
        self.voltages.insert(terminal, voltage);
    }

    /// Set spin accumulation at a terminal
    pub fn set_spin_accumulation(&mut self, terminal: Terminal, mu_s: f64) {
        self.spin_accumulation.insert(terminal, mu_s);
    }

    /// Calculate charge current through a connection
    ///
    /// I = (V1 - V2) / R_total
    pub fn charge_current(&self, connection: &Connection) -> f64 {
        let v1 = self.voltages.get(&connection.from).copied().unwrap_or(0.0);
        let v2 = self.voltages.get(&connection.to).copied().unwrap_or(0.0);
        let r_total = connection.resistor.total_resistance();

        (v1 - v2) / r_total
    }

    /// Calculate spin current through a connection
    ///
    /// I_s = P × I_charge + (μ_s1 - μ_s2) / (2e × R_s)
    pub fn spin_current(&self, connection: &Connection) -> f64 {
        let i_charge = self.charge_current(connection);
        let p = connection.resistor.current_polarization();

        let mu_s1 = self
            .spin_accumulation
            .get(&connection.from)
            .copied()
            .unwrap_or(0.0);
        let mu_s2 = self
            .spin_accumulation
            .get(&connection.to)
            .copied()
            .unwrap_or(0.0);

        // Spin current from charge current polarization
        let i_s_pol = p * i_charge;

        // Spin current from chemical potential difference
        let r_s = (connection.resistor.r_up - connection.resistor.r_down).abs() / 2.0;
        let i_s_diff = if r_s > 0.0 {
            (mu_s1 - mu_s2) / (2.0 * r_s)
        } else {
            0.0
        };

        i_s_pol + i_s_diff
    }

    /// Get all terminals in the circuit
    pub fn terminals(&self) -> Vec<Terminal> {
        let mut terminals = std::collections::HashSet::new();
        for conn in &self.connections {
            terminals.insert(conn.from);
            terminals.insert(conn.to);
        }
        terminals.into_iter().collect()
    }

    /// Calculate total current into a terminal (charge conservation)
    pub fn current_into(&self, terminal: Terminal) -> f64 {
        let mut total = 0.0;

        for conn in &self.connections {
            if conn.to == terminal {
                total += self.charge_current(conn);
            } else if conn.from == terminal {
                total -= self.charge_current(conn);
            }
        }

        total
    }

    /// Calculate total spin current into a terminal
    pub fn spin_current_into(&self, terminal: Terminal) -> f64 {
        let mut total = 0.0;

        for conn in &self.connections {
            if conn.to == terminal {
                total += self.spin_current(conn);
            } else if conn.from == terminal {
                total -= self.spin_current(conn);
            }
        }

        total
    }
}

impl Default for SpinCircuit {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_circuit_creation() {
        let circuit = SpinCircuit::new();
        assert_eq!(circuit.connections.len(), 0);
    }

    #[test]
    fn test_add_connection() {
        let mut circuit = SpinCircuit::new();
        let resistor = SpinResistor::new(100.0, 200.0, 1.0e-6, 1.0e-12);

        circuit.add_connection(Terminal(0), Terminal(1), resistor);
        assert_eq!(circuit.connections.len(), 1);
    }

    #[test]
    fn test_charge_current() {
        let mut circuit = SpinCircuit::new();
        let resistor = SpinResistor::new(100.0, 200.0, 1.0e-6, 1.0e-12);

        circuit.add_connection(Terminal(0), Terminal(1), resistor);
        circuit.set_voltage(Terminal(0), 1.0);
        circuit.set_voltage(Terminal(1), 0.0);

        let i = circuit.charge_current(&circuit.connections[0]);

        // I = V/R, R_total ≈ 66.67 Ω
        assert!(i > 0.01);
        assert!(i < 0.02);
    }

    #[test]
    fn test_terminals() {
        let mut circuit = SpinCircuit::new();
        let resistor = SpinResistor::new(100.0, 200.0, 1.0e-6, 1.0e-12);

        circuit.add_connection(Terminal(0), Terminal(1), resistor.clone());
        circuit.add_connection(Terminal(1), Terminal(2), resistor);

        let terminals = circuit.terminals();
        assert_eq!(terminals.len(), 3);
    }
}
