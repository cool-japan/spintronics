//! Altermagnetic spintronics
//!
//! Altermagnets are a newly classified magnetic phase (2022-2024) that is
//! distinct from both ferromagnets and conventional antiferromagnets. They
//! combine properties previously thought to be mutually exclusive:
//!
//! - **Zero net magnetization** (like antiferromagnets)
//! - **Non-relativistic spin splitting** in k-space (like ferromagnets)
//! - **No spin-orbit coupling required** for spin-split bands
//!
//! ## Symmetry Classification
//!
//! The spin splitting in altermagnets has a characteristic angular symmetry
//! in momentum space, classified by analogy with superconducting order
//! parameters:
//!
//! | Symmetry | Harmonic | Nodes | Example materials |
//! |----------|----------|-------|-------------------|
//! | d-wave   | cos(2φ)  | 4     | RuO2, MnTe, Fe2O3 |
//! | g-wave   | cos(4φ)  | 8     | CrSb              |
//! | i-wave   | cos(6φ)  | 12    | (theoretical)      |
//!
//! ## Transport Properties
//!
//! The non-relativistic spin splitting enables unique transport phenomena:
//!
//! - **Spin-splitter effect**: Spin-dependent conductivity without SOC
//! - **Crystal Hall effect**: Anomalous Hall effect from Berry curvature
//!   of spin-split bands, without spin-orbit coupling
//! - **Spin-charge conversion**: Efficient generation of spin currents
//!
//! ## Physical Background
//!
//! In conventional antiferromagnets, the two sublattices are related by
//! translation symmetry combined with time reversal, which ensures that
//! spin-up and spin-down bands are degenerate everywhere in k-space.
//!
//! In altermagnets, the sublattices are instead related by rotational
//! symmetry combined with time reversal. This weaker constraint allows
//! spin splitting at general k-points while maintaining zero net
//! magnetization. The splitting vanishes only along high-symmetry lines
//! (nodes) determined by the crystal symmetry.
//!
//! ## References
//!
//! - L. Smejkal, J. Sinova, T. Jungwirth, "Emerging Research Landscape of
//!   Altermagnetism", Phys. Rev. X 12, 040501 (2022)
//! - L. Smejkal, J. Sinova, T. Jungwirth, "Beyond Conventional Ferromagnetism
//!   and Antiferromagnetism: A New Phase of Magnetically Ordered Materials",
//!   Phys. Rev. X 12, 031042 (2022)
//! - J. Krempasky et al., "Altermagnetic lifting of Kramers spin degeneracy",
//!   Nature 626, 517-522 (2024)
//!
//! ## Example
//!
//! ```rust
//! use spintronics::altermagnet::{Altermagnet, AltermagneticSymmetry, AltermagnetTransport};
//!
//! // Create a RuO2 altermagnet
//! let ruo2 = Altermagnet::ruo2();
//! assert_eq!(ruo2.symmetry, AltermagneticSymmetry::DWave);
//!
//! // Verify zero net magnetization
//! assert!((ruo2.net_magnetization()).abs() < 1e-10);
//!
//! // Check spin splitting at different angles
//! let splitting_0 = ruo2.splitting_at_angle(0.0);
//! let splitting_45 = ruo2.splitting_at_angle(std::f64::consts::PI / 4.0);
//! assert!(splitting_0.abs() > 1.0); // Maximum at 0 degrees
//! assert!(splitting_45.abs() < 1e-10); // Node at 45 degrees
//!
//! // Compute transport properties
//! let transport = AltermagnetTransport::new(&ruo2, 1.0e6, 1.0e-14)
//!     .expect("Should create transport");
//! let j_s = transport.spin_splitter_current(1.0e10)
//!     .expect("Should compute spin current");
//! assert!(j_s.abs() > 0.0); // Nonzero spin current despite zero magnetization
//! ```

pub mod materials;
pub mod transport;

pub use materials::{Altermagnet, AltermagneticSymmetry};
pub use transport::AltermagnetTransport;
