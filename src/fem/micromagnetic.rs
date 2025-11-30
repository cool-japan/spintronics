//! Micromagnetic FEM solver
//!
//! Solves the Landau-Lifshitz-Gilbert equation using finite elements,
//! particularly useful for complex geometries and material distributions.

use crate::fem::assembly::{assemble_mass_matrix, assemble_stiffness_matrix};
use crate::fem::mesh::Mesh2D;
use crate::fem::solver::{solve_linear_system, SolverType};
use crate::material::Ferromagnet;
use crate::vector3::Vector3;

/// Micromagnetic FEM solver
pub struct MicromagneticFEM {
    /// Finite element mesh
    pub mesh: Mesh2D,

    /// Material properties
    pub material: Ferromagnet,

    /// Magnetization field at each node \[A/m\]
    pub magnetization: Vec<Vector3<f64>>,
}

impl MicromagneticFEM {
    /// Create new micromagnetic FEM solver
    ///
    /// # Arguments
    /// * `mesh` - Finite element mesh
    /// * `material` - Magnetic material
    pub fn new(mesh: Mesh2D, material: Ferromagnet) -> Self {
        let n = mesh.n_nodes();

        // Initialize with saturated magnetization in z-direction
        let m0 = Vector3::new(0.0, 0.0, material.ms);
        let magnetization = vec![m0; n];

        Self {
            mesh,
            material,
            magnetization,
        }
    }

    /// Solve for equilibrium magnetization configuration
    ///
    /// Minimizes the total energy E = E_ex + E_anis + E_zeeman
    ///
    /// # Arguments
    /// * `h_ext` - External magnetic field \[A/m\]
    ///
    /// # Returns
    /// Equilibrium magnetization field
    pub fn solve_equilibrium(&mut self, h_ext: Vector3<f64>) -> Result<Vec<Vector3<f64>>, String> {
        // Assemble system matrices
        let k_stiff = assemble_stiffness_matrix(&self.mesh);
        let _m_mass = assemble_mass_matrix(&self.mesh);

        // Solve for each component separately
        // (Simplified - full micromagnetics requires nonlinear solver)
        let n = self.mesh.n_nodes();

        // For demonstration: solve simplified system
        let m_x = vec![0.0; n];
        let m_y = vec![0.0; n];

        // Apply external field effect (simplified)
        let rhs_z: Vec<f64> = (0..n)
            .map(|_| h_ext.z * self.material.ms / self.material.exchange_a)
            .collect();

        // Solve K * m_z = rhs_z
        let m_z = solve_linear_system(&k_stiff, &rhs_z, SolverType::Jacobi)?;

        // Normalize magnetization to M_s
        for i in 0..n {
            let m = Vector3::new(m_x[i], m_y[i], m_z[i]);
            let mag = m.magnitude();
            if mag > 0.0 {
                self.magnetization[i] = Vector3::new(
                    m.x * self.material.ms / mag,
                    m.y * self.material.ms / mag,
                    m.z * self.material.ms / mag,
                );
            }
        }

        Ok(self.magnetization.clone())
    }

    /// Compute total exchange energy
    ///
    /// E_ex = A_ex ∫ |∇m|² dΩ
    pub fn exchange_energy(&self) -> f64 {
        let k_stiff = assemble_stiffness_matrix(&self.mesh);

        let mut energy = 0.0;

        // E_ex = A_ex * m^T K m (for each component)
        for component in 0..3 {
            let m_comp: Vec<f64> = self
                .magnetization
                .iter()
                .map(|m| match component {
                    0 => m.x,
                    1 => m.y,
                    _ => m.z,
                })
                .collect();

            let km = k_stiff.matvec(&m_comp);
            let comp_energy: f64 = m_comp.iter().zip(km.iter()).map(|(m, k)| m * k).sum();
            energy += comp_energy;
        }

        0.5 * self.material.exchange_a * energy
    }

    /// Compute uniaxial anisotropy energy
    ///
    /// E_anis = -K_u ∫ (m · e_anis)² dΩ
    ///
    /// where e_anis is the easy axis direction
    ///
    /// # Arguments
    /// * `easy_axis` - Uniaxial anisotropy easy axis (normalized)
    pub fn uniaxial_anisotropy_energy(&self, easy_axis: Vector3<f64>) -> f64 {
        let m_mass = assemble_mass_matrix(&self.mesh);

        // Compute (m · e_anis)² for each node
        let m_proj_sq: Vec<f64> = self
            .magnetization
            .iter()
            .map(|m| {
                let proj = m.dot(&easy_axis) / self.material.ms;
                proj * proj
            })
            .collect();

        // Integrate using mass matrix
        let m_proj_sq_integrated = m_mass.matvec(&m_proj_sq);
        let integral: f64 = m_proj_sq
            .iter()
            .zip(m_proj_sq_integrated.iter())
            .map(|(a, b)| a * b)
            .sum();

        -self.material.anisotropy_k * integral
    }

    /// Compute cubic anisotropy energy
    ///
    /// E_cubic = K₁ ∫ (m_x² m_y² + m_y² m_z² + m_z² m_x²) dΩ
    ///
    /// # Arguments
    /// * `k1` - First cubic anisotropy constant \[J/m³\]
    pub fn cubic_anisotropy_energy(&self, k1: f64) -> f64 {
        let m_mass = assemble_mass_matrix(&self.mesh);

        // Compute cubic term for each node
        let cubic_term: Vec<f64> = self
            .magnetization
            .iter()
            .map(|m| {
                let mx2 = (m.x / self.material.ms).powi(2);
                let my2 = (m.y / self.material.ms).powi(2);
                let mz2 = (m.z / self.material.ms).powi(2);
                mx2 * my2 + my2 * mz2 + mz2 * mx2
            })
            .collect();

        // Integrate using mass matrix
        let cubic_integrated = m_mass.matvec(&cubic_term);
        let integral: f64 = cubic_term
            .iter()
            .zip(cubic_integrated.iter())
            .map(|(a, b)| a * b)
            .sum();

        k1 * integral
    }

    /// Compute Zeeman energy
    ///
    /// E_zeeman = -μ₀ ∫ M · H_ext dΩ
    pub fn zeeman_energy(&self, h_ext: Vector3<f64>) -> f64 {
        let m_mass = assemble_mass_matrix(&self.mesh);

        let mu0 = 4.0 * std::f64::consts::PI * 1e-7; // Permeability of free space

        // Compute M · H for each node
        let mh: Vec<f64> = self.magnetization.iter().map(|m| m.dot(&h_ext)).collect();

        // Integrate using mass matrix
        let m_mh = m_mass.matvec(&mh);
        -mu0 * mh.iter().zip(m_mh.iter()).map(|(a, b)| a * b).sum::<f64>()
    }

    /// Total micromagnetic energy (exchange + Zeeman)
    ///
    /// For energy including anisotropy, use `total_energy_with_anisotropy`
    pub fn total_energy(&self, h_ext: Vector3<f64>) -> f64 {
        self.exchange_energy() + self.zeeman_energy(h_ext)
    }

    /// Total micromagnetic energy including uniaxial anisotropy
    ///
    /// # Arguments
    /// * `h_ext` - External field \[A/m\]
    /// * `easy_axis` - Uniaxial anisotropy easy axis (normalized)
    pub fn total_energy_with_anisotropy(
        &self,
        h_ext: Vector3<f64>,
        easy_axis: Vector3<f64>,
    ) -> f64 {
        self.exchange_energy()
            + self.zeeman_energy(h_ext)
            + self.uniaxial_anisotropy_energy(easy_axis)
    }

    /// Complete micromagnetic energy including all terms
    ///
    /// # Arguments
    /// * `h_ext` - External field \[A/m\]
    /// * `easy_axis` - Uniaxial anisotropy easy axis (normalized)
    /// * `demag_factors` - Demagnetization factors (Nx, Ny, Nz)
    pub fn total_energy_complete(
        &self,
        h_ext: Vector3<f64>,
        easy_axis: Vector3<f64>,
        demag_factors: Vector3<f64>,
    ) -> f64 {
        self.exchange_energy()
            + self.zeeman_energy(h_ext)
            + self.uniaxial_anisotropy_energy(easy_axis)
            + self.demagnetization_energy(demag_factors)
    }

    /// Compute demagnetization energy
    ///
    /// E_demag = (μ₀/2) ∫ H_demag · M dΩ
    ///
    /// Note: This is a simplified calculation. Full demagnetization requires
    /// solving Poisson's equation ∇²φ = ∇·M with appropriate boundary conditions.
    /// This implementation provides an approximate demagnetization energy based
    /// on shape anisotropy.
    ///
    /// # Arguments
    /// * `demagnetization_factors` - Nx, Ny, Nz demagnetization factors
    pub fn demagnetization_energy(&self, demagnetization_factors: Vector3<f64>) -> f64 {
        let m_mass = assemble_mass_matrix(&self.mesh);
        let mu0 = 4.0 * std::f64::consts::PI * 1e-7;

        let mut energy = 0.0;

        // Demagnetization field H_d = -N·M
        for component in 0..3 {
            let n_factor = match component {
                0 => demagnetization_factors.x,
                1 => demagnetization_factors.y,
                _ => demagnetization_factors.z,
            };

            let m_comp: Vec<f64> = self
                .magnetization
                .iter()
                .map(|m| match component {
                    0 => m.x,
                    1 => m.y,
                    _ => m.z,
                })
                .collect();

            // E_demag = (μ₀/2) * N * ∫ M² dΩ
            let m_integrated = m_mass.matvec(&m_comp);
            let comp_energy: f64 = m_comp
                .iter()
                .zip(m_integrated.iter())
                .map(|(m, mi)| m * mi)
                .sum();

            energy += n_factor * comp_energy;
        }

        0.5 * mu0 * energy
    }

    /// Compute effective field H_eff from energy derivatives
    ///
    /// H_eff = -(1/μ₀) ∂E/∂M
    ///
    /// # Arguments
    /// * `h_ext` - External applied field \[A/m\]
    /// * `easy_axis` - Anisotropy easy axis (normalized)
    /// * `demag_factors` - Demagnetization factors (Nx, Ny, Nz)
    ///
    /// # Returns
    /// Effective field at each node \[A/m\]
    pub fn compute_effective_field(
        &self,
        h_ext: Vector3<f64>,
        easy_axis: Vector3<f64>,
        demag_factors: Vector3<f64>,
    ) -> Vec<Vector3<f64>> {
        let n = self.mesh.n_nodes();
        let mut h_eff = vec![Vector3::new(0.0, 0.0, 0.0); n];

        // Exchange field: H_ex = (2A_ex/μ₀M_s) ∇²m
        let k_stiff = assemble_stiffness_matrix(&self.mesh);
        let mu0 = 4.0 * std::f64::consts::PI * 1e-7;
        let exchange_coeff = 2.0 * self.material.exchange_a / (mu0 * self.material.ms);

        let inv_ms = 1.0 / self.material.ms;
        for component in 0..3 {
            let m_comp: Vec<f64> = self
                .magnetization
                .iter()
                .map(|m| match component {
                    0 => m.x * inv_ms,
                    1 => m.y * inv_ms,
                    _ => m.z * inv_ms,
                })
                .collect();

            let laplacian = k_stiff.matvec(&m_comp);

            for i in 0..n {
                let h_ex = -exchange_coeff * laplacian[i];
                match component {
                    0 => h_eff[i].x += h_ex,
                    1 => h_eff[i].y += h_ex,
                    _ => h_eff[i].z += h_ex,
                }
            }
        }

        // Anisotropy field: H_anis = (2K_u/μ₀M_s²)(m·e_anis)e_anis
        let anis_coeff =
            2.0 * self.material.anisotropy_k / (mu0 * self.material.ms * self.material.ms);
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            let m_norm = self.magnetization[i] * inv_ms;
            let m_dot_e = m_norm.dot(&easy_axis);
            let h_anis = easy_axis * (anis_coeff * m_dot_e);
            h_eff[i] = h_eff[i] + h_anis;
        }

        // Demagnetization field: H_demag = -N·M
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            let h_demag = Vector3::new(
                -demag_factors.x * self.magnetization[i].x,
                -demag_factors.y * self.magnetization[i].y,
                -demag_factors.z * self.magnetization[i].z,
            );
            h_eff[i] = h_eff[i] + h_demag;
        }

        // External field
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            h_eff[i] = h_eff[i] + h_ext;
        }

        h_eff
    }

    /// Perform one time step of LLG dynamics using semi-implicit scheme
    ///
    /// Solves: dm/dt = -γ₀/(1+α²) [(m × H_eff) + α m × (m × H_eff)]
    ///
    /// # Arguments
    /// * `dt` - Time step \[s\]
    /// * `h_ext` - External field \[A/m\]
    /// * `easy_axis` - Anisotropy easy axis
    /// * `demag_factors` - Demagnetization factors
    ///
    /// # Returns
    /// Maximum change in magnetization (for convergence checking)
    pub fn step_llg(
        &mut self,
        dt: f64,
        h_ext: Vector3<f64>,
        easy_axis: Vector3<f64>,
        demag_factors: Vector3<f64>,
    ) -> f64 {
        let gamma0 = 2.21e5; // Gyromagnetic ratio [m/(A·s)]
        let alpha = self.material.alpha;
        let prefactor = gamma0 / (1.0 + alpha * alpha);

        // Compute effective field
        let h_eff = self.compute_effective_field(h_ext, easy_axis, demag_factors);

        let mut max_change = 0.0;
        let inv_ms = 1.0 / self.material.ms;

        // Update magnetization at each node
        #[allow(clippy::needless_range_loop)]
        for i in 0..self.magnetization.len() {
            let m = self.magnetization[i] * inv_ms; // Normalized
            let h = h_eff[i];

            // LLG terms
            let m_cross_h = m.cross(&h);
            let m_cross_m_cross_h = m.cross(&m_cross_h);

            // dm/dt
            let dmdt = (m_cross_h + m_cross_m_cross_h * alpha) * prefactor;

            // Semi-implicit update: m_new = m + dt * dmdt
            let m_new_unnormalized = m + dmdt * dt;

            // Renormalize to maintain |m| = 1
            let m_new_magnitude = m_new_unnormalized.magnitude();
            let m_new = if m_new_magnitude > 1e-10 {
                m_new_unnormalized * (1.0 / m_new_magnitude)
            } else {
                m
            };

            // Track maximum change
            let change = (m_new - m).magnitude();
            if change > max_change {
                max_change = change;
            }

            // Update magnetization (denormalize back to M_s)
            self.magnetization[i] = m_new * self.material.ms;
        }

        max_change
    }

    /// Run LLG dynamics for a specified time
    ///
    /// # Arguments
    /// * `total_time` - Total simulation time \[s\]
    /// * `dt` - Time step \[s\]
    /// * `h_ext` - External field \[A/m\]
    /// * `easy_axis` - Anisotropy easy axis
    /// * `demag_factors` - Demagnetization factors
    ///
    /// # Returns
    /// Number of steps taken
    pub fn run_dynamics(
        &mut self,
        total_time: f64,
        dt: f64,
        h_ext: Vector3<f64>,
        easy_axis: Vector3<f64>,
        demag_factors: Vector3<f64>,
    ) -> usize {
        let n_steps = (total_time / dt) as usize;

        for _step in 0..n_steps {
            self.step_llg(dt, h_ext, easy_axis, demag_factors);
        }

        n_steps
    }

    /// Get average magnetization
    pub fn average_magnetization(&self) -> Vector3<f64> {
        let n = self.magnetization.len() as f64;
        if n == 0.0 {
            return Vector3::new(0.0, 0.0, 0.0);
        }

        let mut sum = Vector3::new(0.0, 0.0, 0.0);
        for m in &self.magnetization {
            sum = sum + *m;
        }

        Vector3::new(sum.x / n, sum.y / n, sum.z / n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_micromagnetic_fem_creation() {
        let mesh = Mesh2D::rectangle(100e-9, 50e-9, 10e-9).unwrap();
        let material = Ferromagnet::permalloy();

        let fem = MicromagneticFEM::new(mesh, material);

        assert_eq!(fem.magnetization.len(), fem.mesh.n_nodes());
        assert!(fem.magnetization[0].magnitude() > 0.0);
    }

    #[test]
    fn test_average_magnetization() {
        let mesh = Mesh2D::rectangle(100e-9, 50e-9, 20e-9).unwrap();
        let material = Ferromagnet::permalloy();
        let ms = material.ms;

        let fem = MicromagneticFEM::new(mesh, material);
        let m_avg = fem.average_magnetization();

        // Should be close to (0, 0, M_s)
        assert!(m_avg.z > 0.0);
        assert!((m_avg.magnitude() - ms).abs() < ms * 0.1);
    }

    #[test]
    fn test_energy_methods_exist() {
        let mesh = Mesh2D::rectangle(50e-9, 50e-9, 10e-9).unwrap();
        let material = Ferromagnet::yig();

        let fem = MicromagneticFEM::new(mesh, material);
        let h_ext = Vector3::new(0.0, 0.0, 1000.0);

        // Test that energy methods can be called (values may need refinement)
        let _e_ex = fem.exchange_energy();
        let _e_z = fem.zeeman_energy(h_ext);
        let _e_total = fem.total_energy(h_ext);

        // Just verify the methods exist and don't panic on small mesh
    }

    #[test]
    fn test_uniaxial_anisotropy_energy() {
        let mesh = Mesh2D::rectangle(50e-9, 50e-9, 10e-9).unwrap();
        let material = Ferromagnet::permalloy();

        let fem = MicromagneticFEM::new(mesh, material);

        // Easy axis along z
        let easy_axis = Vector3::new(0.0, 0.0, 1.0);

        // Calculate anisotropy energy
        let e_anis = fem.uniaxial_anisotropy_energy(easy_axis);

        // Should be negative (energy is minimized along easy axis)
        assert!(
            e_anis <= 0.0,
            "Uniaxial anisotropy energy should be non-positive"
        );
    }

    #[test]
    fn test_cubic_anisotropy_energy() {
        let mesh = Mesh2D::rectangle(50e-9, 50e-9, 10e-9).unwrap();
        let material = Ferromagnet::permalloy();

        let fem = MicromagneticFEM::new(mesh, material);

        // First cubic anisotropy constant
        let k1 = -5e3; // J/m³ (typical for Fe)

        // Calculate cubic anisotropy energy
        let _e_cubic = fem.cubic_anisotropy_energy(k1);

        // Just verify it doesn't panic
    }

    #[test]
    fn test_total_energy_with_anisotropy() {
        let mesh = Mesh2D::rectangle(50e-9, 50e-9, 15e-9).unwrap();
        let material = Ferromagnet::permalloy();

        let fem = MicromagneticFEM::new(mesh, material);

        let h_ext = Vector3::new(0.0, 0.0, 1000.0);
        let easy_axis = Vector3::new(0.0, 0.0, 1.0);

        // Calculate total energy with anisotropy
        let _e_total_anis = fem.total_energy_with_anisotropy(h_ext, easy_axis);

        // Just verify the method doesn't panic
    }

    #[test]
    fn test_demagnetization_energy() {
        let mesh = Mesh2D::rectangle(100e-9, 50e-9, 10e-9).unwrap();
        let material = Ferromagnet::permalloy();

        let fem = MicromagneticFEM::new(mesh, material);

        // Thin film demagnetization factors (approximate)
        // For a thin film: Nx ≈ 0, Ny ≈ 0, Nz ≈ 1
        let demag_factors = Vector3::new(0.0, 0.0, 1.0);

        // Calculate demagnetization energy
        let e_demag = fem.demagnetization_energy(demag_factors);

        // Should be positive (demagnetization costs energy)
        assert!(
            e_demag >= 0.0,
            "Demagnetization energy should be non-negative"
        );
    }

    #[test]
    fn test_total_energy_complete() {
        let mesh = Mesh2D::rectangle(50e-9, 50e-9, 15e-9).unwrap();
        let material = Ferromagnet::permalloy();

        let fem = MicromagneticFEM::new(mesh, material);

        let h_ext = Vector3::new(0.0, 0.0, 1000.0);
        let easy_axis = Vector3::new(0.0, 0.0, 1.0);
        let demag_factors = Vector3::new(0.0, 0.0, 1.0);

        // Calculate complete energy
        let _e_complete = fem.total_energy_complete(h_ext, easy_axis, demag_factors);

        // Just verify the method doesn't panic
    }

    #[test]
    fn test_effective_field_calculation() {
        let mesh = Mesh2D::rectangle(50e-9, 50e-9, 10e-9).unwrap();
        let material = Ferromagnet::permalloy();

        let fem = MicromagneticFEM::new(mesh, material);

        let h_ext = Vector3::new(1000.0, 0.0, 0.0);
        let easy_axis = Vector3::new(0.0, 0.0, 1.0);
        let demag_factors = Vector3::new(0.0, 0.0, 1.0);

        let h_eff = fem.compute_effective_field(h_ext, easy_axis, demag_factors);

        // Should have field at each node
        assert_eq!(h_eff.len(), fem.mesh.n_nodes());

        // Fields should be non-zero
        assert!(h_eff.iter().any(|h| h.magnitude() > 0.0));
    }

    #[test]
    fn test_llg_single_step() {
        let mesh = Mesh2D::rectangle(50e-9, 50e-9, 10e-9).unwrap();
        let material = Ferromagnet::permalloy();

        let mut fem = MicromagneticFEM::new(mesh, material);

        let dt = 1e-13; // 0.1 ps
        let h_ext = Vector3::new(10000.0, 0.0, 0.0);
        let easy_axis = Vector3::new(0.0, 0.0, 1.0);
        let demag_factors = Vector3::new(0.0, 0.0, 1.0);

        let m_initial = fem.average_magnetization();

        // Take one LLG step
        let max_change = fem.step_llg(dt, h_ext, easy_axis, demag_factors);

        let m_after = fem.average_magnetization();

        // Magnetization should have changed
        assert!(max_change > 0.0);

        // Magnetization magnitude should be conserved
        let mag_initial = m_initial.magnitude();
        let mag_after = m_after.magnitude();
        assert!((mag_initial - mag_after).abs() / mag_initial < 0.01);
    }

    #[test]
    fn test_llg_dynamics_run() {
        let mesh = Mesh2D::rectangle(50e-9, 50e-9, 10e-9).unwrap();
        let material = Ferromagnet::permalloy();
        let ms = material.ms; // Extract before move

        let mut fem = MicromagneticFEM::new(mesh, material);

        let total_time = 1e-12; // 1 ps
        let dt = 1e-13; // 0.1 ps
        let h_ext = Vector3::new(5000.0, 0.0, 0.0);
        let easy_axis = Vector3::new(0.0, 0.0, 1.0);
        let demag_factors = Vector3::new(0.0, 0.0, 1.0);

        let n_steps = fem.run_dynamics(total_time, dt, h_ext, easy_axis, demag_factors);

        // Should have taken 10 steps
        assert_eq!(n_steps, 10);

        // Magnetization should still be normalized
        let m_avg = fem.average_magnetization();
        assert!((m_avg.magnitude() - ms).abs() / ms < 0.1);
    }

    #[test]
    fn test_llg_dynamics_stability() {
        // Test that LLG dynamics doesn't blow up
        let mesh = Mesh2D::rectangle(50e-9, 50e-9, 10e-9).unwrap();
        let material = Ferromagnet::permalloy();

        let mut fem = MicromagneticFEM::new(mesh, material);

        let h_ext = Vector3::new(5000.0, 0.0, 0.0);
        let easy_axis = Vector3::new(0.0, 0.0, 1.0);
        let demag_factors = Vector3::new(0.0, 0.0, 0.0);

        // Run a few steps
        for _ in 0..10 {
            let _change = fem.step_llg(1e-13, h_ext, easy_axis, demag_factors);
        }

        // Verify magnetization is still reasonable
        let m_avg = fem.average_magnetization();
        assert!(m_avg.magnitude() > 0.0);
        assert!(m_avg.magnitude().is_finite());
    }
}
