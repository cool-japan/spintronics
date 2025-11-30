//! Experimental Benchmarks and Validation
//!
//! This module provides reference values and benchmark functions to validate
//! simulation results against known experimental data from the literature.
//!
//! Each benchmark includes:
//! - Experimental reference values
//! - Material parameters
//! - Measurement conditions
//! - Citation to original paper
//!
//! ## Usage
//!
//! Benchmarks can be used to:
//! - Validate your simulation parameters
//! - Compare computed values with experiments
//! - Ensure physical correctness

// Note: imports for future use in more complex benchmarks
#[allow(unused_imports)]
use crate::constants::*;
#[allow(unused_imports)]
use crate::vector3::Vector3;

/// Experimental benchmark data
#[derive(Debug, Clone)]
pub struct Benchmark {
    /// Name of the experiment/effect
    pub name: String,

    /// Reference citation
    pub citation: String,

    /// Experimental value
    pub experimental_value: f64,

    /// Unit of the value
    pub unit: String,

    /// Experimental uncertainty (if known)
    pub uncertainty: Option<f64>,

    /// Temperature at which measurement was taken \[K\]
    pub temperature: f64,

    /// Additional notes
    pub notes: String,
}

impl Benchmark {
    /// Compare computed value with experimental benchmark
    ///
    /// Returns relative error: (computed - experimental) / experimental
    pub fn compare(&self, computed: f64) -> f64 {
        (computed - self.experimental_value) / self.experimental_value
    }

    /// Check if computed value agrees within experimental uncertainty
    pub fn agrees(&self, computed: f64) -> bool {
        if let Some(unc) = self.uncertainty {
            (computed - self.experimental_value).abs() <= unc
        } else {
            // If no uncertainty given, use 10% as reasonable agreement
            self.compare(computed).abs() < 0.1
        }
    }
}

/// Saitoh et al. 2006: First observation of spin pumping and ISHE
///
/// "Conversion of spin current into charge current at room temperature:
/// Inverse spin-Hall effect", Applied Physics Letters 88, 182509 (2006)
pub fn saitoh_2006_ishe() -> Benchmark {
    Benchmark {
        name: "Spin pumping voltage (Py/Pt)".to_string(),
        citation: "Saitoh et al., APL 88, 182509 (2006)".to_string(),
        experimental_value: 0.72e-6, // 0.72 µV
        unit: "V".to_string(),
        uncertainty: Some(0.05e-6),
        temperature: 300.0,
        notes: "FMR at 9.42 GHz, Py(20nm)/Pt(10nm), RF power 200 mW".to_string(),
    }
}

/// Uchida et al. 2008: First observation of spin Seebeck effect
///
/// "Observation of the spin Seebeck effect", Nature 455, 778 (2008)
pub fn uchida_2008_sse() -> Benchmark {
    Benchmark {
        name: "Spin Seebeck voltage (LaY₂Fe₅O₁₂/Pt)".to_string(),
        citation: "Uchida et al., Nature 455, 778 (2008)".to_string(),
        experimental_value: 10e-6, // ~10 µV for 10K temperature difference
        unit: "V".to_string(),
        uncertainty: Some(2e-6),
        temperature: 300.0,
        notes: "ΔT = 10 K across LaY₂Fe₅O₁₂ film with Pt detector".to_string(),
    }
}

/// Liu et al. 2012: Spin-torque switching by spin Hall effect
///
/// "Spin-Torque Switching with the Giant Spin Hall Effect of Tantalum",
/// Science 336, 555 (2012)
pub fn liu_2012_sot() -> Benchmark {
    Benchmark {
        name: "Critical switching current density (Ta/CoFeB)".to_string(),
        citation: "Liu et al., Science 336, 555 (2012)".to_string(),
        experimental_value: 2.2e11, // 2.2 × 10¹¹ A/m²
        unit: "A/m²".to_string(),
        uncertainty: Some(0.3e11),
        temperature: 300.0,
        notes: "Ta(3nm)/CoFeB(0.9nm)/MgO, spin Hall angle θ_SH ≈ -0.15".to_string(),
    }
}

/// Nagaosa & Tokura 2013: Skyrmion size in helimagnets
///
/// "Topological properties and dynamics of magnetic skyrmions",
/// Nature Nanotechnology 8, 899 (2013)
pub fn skyrmion_size_fege() -> Benchmark {
    Benchmark {
        name: "Skyrmion diameter (FeGe)".to_string(),
        citation: "Nagaosa & Tokura, Nat. Nanotech. 8, 899 (2013)".to_string(),
        experimental_value: 70e-9, // 70 nm
        unit: "m".to_string(),
        uncertainty: Some(10e-9),
        temperature: 260.0, // Below T_C = 278 K
        notes: "B3-type helimagnet, measured by Lorentz TEM".to_string(),
    }
}

/// Johnson & Silsbee 1985: Spin diffusion length in metals
///
/// "Interfacial charge-spin coupling: Injection and detection of spin
/// magnetization in metals", Physical Review Letters 55, 1790 (1985)
pub fn spin_diffusion_length_al() -> Benchmark {
    Benchmark {
        name: "Spin diffusion length (Aluminum)".to_string(),
        citation: "Johnson & Silsbee, PRL 55, 1790 (1985)".to_string(),
        experimental_value: 600e-9, // 600 nm
        unit: "m".to_string(),
        uncertainty: Some(100e-9),
        temperature: 4.2, // Liquid helium temperature
        notes: "Pure aluminum at low temperature".to_string(),
    }
}

/// Gilbert damping parameter for common materials
pub fn gilbert_damping_py() -> Benchmark {
    Benchmark {
        name: "Gilbert damping (Permalloy)".to_string(),
        citation: "Multiple sources, consensus value".to_string(),
        experimental_value: 0.01,
        unit: "dimensionless".to_string(),
        uncertainty: Some(0.002),
        temperature: 300.0,
        notes: "Ni₈₀Fe₂₀, one of the lowest damping metallic ferromagnets".to_string(),
    }
}

/// Gilbert damping for YIG (ultra-low damping)
pub fn gilbert_damping_yig() -> Benchmark {
    Benchmark {
        name: "Gilbert damping (YIG)".to_string(),
        citation: "Multiple sources, typical value".to_string(),
        experimental_value: 0.0001,
        unit: "dimensionless".to_string(),
        uncertainty: Some(0.00005),
        temperature: 300.0,
        notes: "Y₃Fe₅O₁₂, lowest known damping at room temperature".to_string(),
    }
}

/// Spin Hall angle for Pt
pub fn spin_hall_angle_pt() -> Benchmark {
    Benchmark {
        name: "Spin Hall angle (Platinum)".to_string(),
        citation: "Multiple measurements, typical value".to_string(),
        experimental_value: 0.07,
        unit: "dimensionless".to_string(),
        uncertainty: Some(0.02), // Varies significantly between samples
        temperature: 300.0,
        notes: "Widely used for spin-charge conversion, value depends on purity".to_string(),
    }
}

/// Spin Hall angle for Ta (giant SHE)
pub fn spin_hall_angle_ta() -> Benchmark {
    Benchmark {
        name: "Spin Hall angle (Tantalum β-phase)".to_string(),
        citation: "Liu et al., Science 336, 555 (2012)".to_string(),
        experimental_value: -0.15,
        unit: "dimensionless".to_string(),
        uncertainty: Some(0.03),
        temperature: 300.0,
        notes: "β-Ta phase shows giant spin Hall effect, negative sign".to_string(),
    }
}

/// DMI constant for Pt/CoFeB interface
pub fn dmi_constant_pt_cofeb() -> Benchmark {
    Benchmark {
        name: "DMI constant (Pt/CoFeB interface)".to_string(),
        citation: "Multiple sources, typical value".to_string(),
        experimental_value: 1.2e-3, // 1.2 mJ/m²
        unit: "J/m²".to_string(),
        uncertainty: Some(0.3e-3),
        temperature: 300.0,
        notes: "Interfacial DMI, stabilizes Néel skyrmions".to_string(),
    }
}

/// Collection of all benchmarks
pub fn all_benchmarks() -> Vec<Benchmark> {
    vec![
        saitoh_2006_ishe(),
        uchida_2008_sse(),
        liu_2012_sot(),
        skyrmion_size_fege(),
        spin_diffusion_length_al(),
        gilbert_damping_py(),
        gilbert_damping_yig(),
        spin_hall_angle_pt(),
        spin_hall_angle_ta(),
        dmi_constant_pt_cofeb(),
    ]
}

/// Print all benchmarks in a formatted table
pub fn print_benchmark_table() {
    println!("=== Experimental Benchmarks ===\n");
    println!("{:<40} {:<15} {:<30}", "Quantity", "Value", "Reference");
    println!("{}", "-".repeat(85));

    for bench in all_benchmarks() {
        let value_str = if let Some(unc) = bench.uncertainty {
            format!(
                "{:.2e} ± {:.2e} {}",
                bench.experimental_value, unc, bench.unit
            )
        } else {
            format!("{:.2e} {}", bench.experimental_value, bench.unit)
        };

        println!(
            "{:<40} {:<15} {:<30}",
            bench.name,
            value_str,
            bench.citation.split(',').next().unwrap_or("")
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_benchmark_compare() {
        let bench = gilbert_damping_py();
        let computed = 0.011;

        let rel_error = bench.compare(computed);
        assert!((rel_error - 0.1).abs() < 1e-10); // (0.011 - 0.01)/0.01 = 0.1
    }

    #[test]
    fn test_benchmark_agrees() {
        let bench = gilbert_damping_py();

        // Should agree within uncertainty
        assert!(bench.agrees(0.011));
        assert!(bench.agrees(0.009));

        // Should not agree (too far)
        assert!(!bench.agrees(0.020));
    }

    #[test]
    fn test_all_benchmarks() {
        let benches = all_benchmarks();
        assert!(benches.len() >= 10);

        // Check that all have proper citations
        for bench in benches {
            assert!(!bench.citation.is_empty());
            assert!(bench.experimental_value != 0.0);
        }
    }
}
