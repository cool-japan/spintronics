//! Finite Element Method (FEM) for Micromagnetic Simulations
//!
//! This module provides FEM capabilities for solving micromagnetic problems,
//! particularly useful for complex geometries and material distributions where
//! finite-difference methods become impractical.
//!
//! ## Features
//!
//! - **Mesh Generation**: Delaunay triangulation for 2D/3D domains
//! - **Sparse Linear Solvers**: CG, BiCGSTAB, AMG for efficient solution
//! - **Micromagnetic FEM**: Specialized elements for magnetization dynamics
//! - **Parallel Assembly**: Thread-parallel stiffness matrix assembly
//!
//! ## Physical Background
//!
//! The micromagnetic energy functional is discretized using FEM:
//!
//! ```text
//! E_total = E_exchange + E_anisotropy + E_zeeman + E_demag
//! ```
//!
//! For the LLG equation, we solve:
//! ```text
//! ∂m/∂t = -γ (m × H_eff) + α (m × ∂m/∂t)
//! ```
//!
//! ## Example
//!
//! ```rust,ignore
//! use spintronics::fem::{Mesh2D, MicromagneticFEM};
//! use spintronics::material::Ferromagnet;
//!
//! // Generate triangular mesh
//! let mesh = Mesh2D::rectangle(100e-9, 50e-9, 1e-9)?;
//!
//! // Setup micromagnetic problem
//! let material = Ferromagnet::permalloy();
//! let fem = MicromagneticFEM::new(mesh, material);
//!
//! // Solve for equilibrium magnetization
//! let m_eq = fem.solve_equilibrium()?;
//! ```
//!
//! ## References
//!
//! - A. Abert et al., "A self-consistent spin-diffusion model for micromagnetics",
//!   Scientific Reports 6, 16/2016
//! - A. Vansteenkiste et al., "The design and verification of MuMax3",
//!   AIP Advances 4, 107133 (2014)

pub mod assembly;
pub mod element;
pub mod mesh;
pub mod micromagnetic;
pub mod solver;

pub use assembly::{
    assemble_mass_matrix, assemble_mass_matrix_parallel, assemble_stiffness_matrix,
    assemble_stiffness_matrix_parallel,
};
pub use element::{TetrahedralElement, TriangularElement};
pub use mesh::{Element2D, Mesh2D, Mesh3D, Node};
pub use micromagnetic::MicromagneticFEM;
pub use solver::{
    solve_linear_system, solve_linear_system_with_params, Preconditioner, SolverParams, SolverType,
};
