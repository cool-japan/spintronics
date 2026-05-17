//! # spintronics
//!
//! **Version 0.6.0** - Spin wave extensions, HOTI, axion electrodynamics, data formats, ML autodiff
//!
//! A pure Rust library for simulating spin dynamics, spin current generation,
//! and conversion phenomena in magnetic materials and topological materials.
//!
//! **Copyright (c) 2025 COOLJAPAN OÜ (Team KitaSan)**
//!
//! Licensed under MIT OR Apache-2.0
//!
//! ## Overview
//!
//! This library implements physical models established by Prof. Eiji Saitoh's
//! research group (Univ. Tokyo/RIKEN) and the broader spintronics community:
//!
//! ### Core Physics Effects
//! - **Spin Pumping**: Generation of spin current from magnetization precession
//! - **Inverse Spin Hall Effect (ISHE)**: Conversion of spin current to charge current
//! - **Spin Seebeck Effect (SSE)**: Thermal generation of spin current
//! - **Spin-Orbit Torque (SOT)**: Current-driven magnetization switching
//! - **Dzyaloshinskii-Moriya Interaction (DMI)**: Skyrmion stabilization
//! - **Topological Hall Effect**: Berry phase from skyrmion textures
//! - **Rashba Effect**: 2DEG spin splitting and spin-momentum locking
//! - **Edelstein Effect**: Spin-to-charge conversion in non-centrosymmetric systems
//! - **Spin Nernst Effect**: Thermal gradient → transverse spin current
//!
//! ### Key Features
//! - ✅ **1200 tests passing** (1200 library tests, v0.6.0)
//! - ✅ **Spin wave extensions** - Damon-Eshbach non-reciprocity, BVMSW negative v_g, surface spin waves, spectral magnon solver (v0.6.0)
//! - ✅ **Higher-order topological insulators** - BBH model, breathing kagome, corner states, nested Wilson loops (v0.6.0)
//! - ✅ **Axion electrodynamics** - 3D Berry-curvature θ-term, axion magnon-photon coupling, topological ME polarizability (v0.6.0)
//! - ✅ **Data export** - VTI, XDMF, pure-Rust NetCDF3 Classic, Zarr v2 writers (v0.6.0)
//! - ✅ **ML autodiff** - Reverse-mode AD tape, SGD/Adam/L-BFGS optimizers, differentiable physics functions (v0.6.0)
//! - ✅ **41 examples** organized by difficulty (v0.6.0)
//! - ✅ **Non-collinear magnetism** - Spin spirals (cycloidal/helical/conical), Luttinger-Tisza ground state (v0.5.0)
//! - ✅ **Multiferroics** - Magnetoelectric tensor, KNB mechanism, DM polarization, toroidal moments (v0.5.0)
//! - ✅ **Kane-Mele QSH** - Z2 topological invariant, Rashba/staggered potential, helical edge states (v0.5.0)
//! - ✅ **Nonlinear magnon physics** - Four-magnon scattering, Suhl instability, parametric amplification (v0.5.0)
//! - ✅ **Quantum magnonics** - Holstein-Primakoff, Bogoliubov transform, zero-point fluctuations (v0.4.0)
//! - ✅ **NEGF transport** - Keldysh formalism, Landauer-Büttiker, shot noise, spin accumulation (v0.4.0)
//! - ✅ **Topological magnon bands** - Haldane model, Chern numbers, Berry curvature, edge modes (v0.4.0)
//! - ✅ **Cavity extensions** - Tavis-Cummings, polaritons, optomagnonics, frequency combs (v0.4.0)
//! - ✅ **Random anisotropy** - Imry-Ma disorder, Harris criterion, LLG coupling (v0.4.0)
//! - ✅ **Math primitives** - Complex, CMatrix with Gauss-Jordan inverse + TQLI eigendecomposition (v0.4.0)
//! - ✅ **Advanced integrators** - Dormand-Prince RK5(4)/RK8(7), symplectic methods (v0.3.0)
//! - ✅ **New physics modules** - Altermagnets, Orbitronics, Frustrated magnets, Hopfions, Magnon BEC, Magnetoelastics (v0.3.0)
//! - ✅ **Interactive web demo** - Axum + HTMX subcrate with 4 physics simulations (v0.2.0)
//! - ✅ **5 experimental validations** against landmark papers
//! - ✅ **WebAssembly support** for browser-based simulations
//! - ✅ **FEM solver** with advanced iterative methods
//! - ✅ **Performance optimized** - 21 inline attributes on hot-path functions (v0.2.0)
//! - ✅ **Memory optimized** - Pool allocator reduces allocations by 99% (v0.2.0)
//! - ✅ **Python bindings** via PyO3 (v0.2.0)
//! - ✅ **Serde serialization** for data interchange (v0.2.0)
//! - ✅ **HDF5 export** for large datasets (v0.2.0)
//! - ✅ **Unit validation** - 14 validators for physical quantities (v0.2.0)
//! - ✅ **Zero warnings** - production-quality code
//!
//! ## Key References
//!
//! - E. Saitoh et al., "Conversion of spin current into charge current at room
//!   temperature: Inverse spin-Hall effect", *Appl. Phys. Lett.* **88**, 182509 (2006)
//! - K. Uchida et al., "Observation of the spin Seebeck effect",
//!   *Nature* **455**, 778-781 (2008)
//! - I. M. Miron et al., "Perpendicular switching of a single ferromagnetic layer
//!   induced by in-plane current injection", *Nature* **476**, 189 (2011)
//! - S. Woo et al., "Observation of room-temperature magnetic skyrmions",
//!   *Nat. Mater.* **15**, 501 (2016)
//!
//! ## Architecture
//!
//! The library is organized into ~28 physics-focused modules:
//!
//! ### Core Infrastructure
//! - [`constants`]: Physical constants (ℏ, γ, e, μ_B, k_B, 20+ NIST-validated values)
//! - [`math`]: Math primitives — `Complex`, `CMatrix` (Gauss-Jordan + TQLI eigendecomposition, v0.4.0)
//! - [`vector3`]: Optimized 3D vector operations for spin/magnetization
//! - [`units`]: Unit validation - 14 validators for physical quantities (v0.2.0)
//! - [`error`]: Error handling and result types
//!
//! ### Materials & Properties
//! - [`material`]: Material properties (ferromagnets, interfaces, 2D magnets, topological insulators, Weyl semimetals)
//!
//! ### Dynamics & Transport
//! - [`dynamics`]: Time evolution solvers (LLG equation, RK4, Heun, adaptive methods)
//! - [`transport`]: Spin transport phenomena (spin pumping, diffusion)
//!
//! ### Physical Effects
//! - [`effect`]: Spin-charge conversion (ISHE, SSE, SOT, Rashba, Edelstein, Spin Nernst, Topological Hall)
//! - [`magnon`]: Magnon propagation and spin wave dynamics
//! - [`thermo`]: Thermoelectric effects (ANE, thermal magnon transport, multilayers)
//! - [`texture`]: Magnetic textures (skyrmions, domain walls, DMI, topological charge)
//! - [`spinwave`]: Spin wave theory — dispersion, modes, quantization, Damon-Eshbach, BVMSW, surface waves (v0.6.0)
//!
//! ### Specialized Physics
//! - [`afm`]: Antiferromagnetic dynamics for THz spintronics
//! - [`stochastic`]: Thermal fluctuations and finite-temperature effects
//! - [`cavity`]: Cavity magnonics - Tavis-Cummings, polaritons, optomagnonics (v0.3.0/v0.4.0)
//! - [`altermagnet`]: Altermagnetic materials (RuO2, CrSb, MnTe; spin-splitter effect, v0.3.0)
//! - [`orbitronics`]: Orbital Hall effect, orbital torques (v0.3.0)
//! - [`frustrated`]: Frustrated magnets, spin ice, kagome (v0.3.0)
//! - [`quantum`]: Quantum magnonics — Holstein-Primakoff, Bogoliubov, zero-point fluctuations (v0.4.0)
//! - [`negf`]: Non-equilibrium transport — Keldysh, Landauer-Büttiker, shot noise (v0.4.0)
//! - [`topomagnon`]: Topological magnon bands, HOTI corner states, nested Wilson loops, axion electrodynamics (v0.4.0–v0.6.0)
//! - [`noncollinear`]: Non-collinear magnetism — spin spirals, Luttinger-Tisza method (v0.5.0)
//! - [`multiferroic`]: Magnetoelectric coupling — ME tensor, KNB mechanism, DM polarization (v0.5.0)
//! - [`autodiff`]: Reverse-mode autodiff — tape, SGD/Adam/L-BFGS optimizers, differentiable physics (v0.6.0, feature-gated)
//!
//! ### Coupled Systems
//! - [`circuit`]: Spin circuit elements (resistors, networks, spin accumulation)
//! - [`fluid`]: Spin-vorticity coupling in liquid metals (Barnett effect)
//! - [`mech`]: Nanomechanical spintronics (Barnett, Einstein-de Haas, cantilever coupling)
//! - [`ai`]: Physical reservoir computing with magnon dynamics
//!
//! ### Computational Tools
//! - [`fem`]: Finite element method (Delaunay mesh, iterative solvers, micromagnetics)
//! - [`memory`]: Memory pool allocator for high-performance simulations (v0.2.0)
//! - [`simd`]: SIMD-friendly batch processing for vector operations (v0.3.0)
//! - [`parallel`]: Domain decomposition and parallel sweeps (v0.3.0, feature-gated)
//! - [`builder`]: Type-state SimulationBuilder for validated construction (v0.3.0)
//!
//! ### Data & Validation
//! - [`visualization`]: Data export (VTK, CSV, JSON, HDF5, VTI, XDMF, NetCDF3, Zarr v2, v0.6.0)
//! - [`validation`]: Experimental validation tests against landmark papers
//! - `python`: Python bindings via PyO3 (v0.2.0, optional feature)
//!
//! ## Quick Start
//!
//! ```rust
//! use spintronics::prelude::*;
//!
//! // Setup materials (YIG/Pt system)
//! let yig = Ferromagnet::yig();
//! let interface = SpinInterface::yig_pt();
//! let pt_strip = InverseSpinHall::platinum();
//!
//! // Initialize magnetization state
//! let m = Vector3::new(1.0, 0.0, 0.0);
//! let h_ext = Vector3::new(0.0, 0.0, 1.0);
//!
//! // Solve LLG equation
//! let dm_dt = calc_dm_dt(m, h_ext, GAMMA, yig.alpha);
//!
//! // Calculate spin pumping current
//! let js = spin_pumping_current(&interface, m, dm_dt);
//!
//! // Validate physical quantities (v0.2.0)
//! assert!(is_valid_magnetization(yig.ms));
//! assert!(is_valid_damping(yig.alpha));
//!
//! // Convert to electric field via ISHE
//! let e_field = pt_strip.convert(interface.normal, js);
//! ```

#![warn(missing_docs)]
#![warn(clippy::all)]
#![allow(clippy::doc_lazy_continuation)]
#![allow(clippy::doc_markdown)]

pub mod afm;
pub mod altermagnet;
pub mod caloritronics;
pub mod math;

#[cfg(not(target_arch = "wasm32"))]
pub mod ai;

pub mod benchmark;
pub mod builder;
pub mod cavity;
pub mod circuit;
pub mod constants;
pub mod dynamics;
pub mod effect;
pub mod error;

#[cfg(feature = "fem")]
pub mod fem;

pub mod fluid;
pub mod frustrated;
pub mod io;
pub mod llg;

#[cfg(not(target_arch = "wasm32"))]
pub mod magnon;

pub mod material;
pub mod mech;
pub mod memory;
pub mod multiferroic;

#[cfg(feature = "parallel")]
pub mod parallel;

pub mod orbitronics;
pub mod simd;
pub mod spinwave;

#[cfg(feature = "scirs2")]
pub mod stochastic;

pub mod negf;
pub mod noncollinear;
pub mod quantum;
pub mod texture;
pub mod thermo;
pub mod topomagnon;
pub mod transport;
pub mod units;
pub mod validation;
pub mod vector3;
pub mod visualization;

#[cfg(feature = "wasm")]
pub mod wasm;

#[cfg(feature = "python")]
pub mod python;

#[cfg(feature = "autodiff")]
pub mod autodiff;

pub mod prelude;

pub use vector3::Vector3;
