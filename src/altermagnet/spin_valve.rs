//! Altermagnet spin valve: giant magnetoresistance without ferromagnetism
//!
//! Conventional giant magnetoresistance (GMR) requires two *ferromagnetic*
//! layers whose net moments can be aligned (parallel) or reversed
//! (antiparallel) by an external field, modulating the spin-dependent
//! scattering seen by conduction electrons. This module implements the
//! altermagnetic analog: two altermagnetic layers, each with **zero net
//! magnetization**, whose non-relativistic k-space spin splitting still
//! produces a strong, angle-dependent magnetoresistance.
//!
//! ## Mechanism
//!
//! Each altermagnetic layer spin-splits its bands with an angular pattern
//! set by its `AltermagneticSymmetry` harmonic order `n` (2, 4, or 6 for
//! d-, g-, and i-wave order): the splitting is proportional to
//! `cos(n * phi)`, where `phi` is measured from the layer's Neel/crystal
//! axis. Rotating one layer's crystal axis by `theta` relative to the other
//! rotates its spin-split "lobes" in momentum space by the same angle.
//!
//! - At `theta = 0` the two layers' lobes and nodes are perfectly aligned:
//!   the same momentum-space channels are transmissive in both layers, so
//!   the spin-dependent scattering mismatch -- and hence the
//!   magnetoresistance -- is minimal.
//! - At the **channel-swap angle** `theta = pi/n`, the top layer's lobes are
//!   rotated exactly onto the bottom layer's nodes (and vice versa): the
//!   channels transmissive in one layer are now blocked in the other,
//!   maximizing the scattering mismatch and the resistance.
//!
//! This gives the same functional form as conventional GMR,
//! `(1 - cos(n*theta)) / 2`, with the ordinary ferromagnetic case recovered
//! formally as `n = 1`; altermagnets instead have `n in {2, 4, 6}` and,
//! crucially, **no ferromagnetic moment ever appears** -- both layers stay
//! antiferromagnetically ordered (equal, opposite sublattice moments) at
//! every relative angle.
//!
//! ## References
//!
//! - L. Smejkal et al., "Emerging Research Landscape of Altermagnetism",
//!   Phys. Rev. X 12, 040501 (2022)
//! - L. Smejkal et al., "Beyond Conventional Ferromagnetism and
//!   Antiferromagnetism", Phys. Rev. X 12, 031042 (2022)
//! - M. N. Baibich et al., "Giant Magnetoresistance of (001)Fe/(001)Cr
//!   Magnetic Superlattices", Phys. Rev. Lett. 61, 2472 (1988) (conventional
//!   GMR, the effect generalized here)
//!
//! ## Example
//!
//! ```rust
//! use spintronics::altermagnet::AltermagnetSpinValve;
//! use std::f64::consts::PI;
//!
//! let valve = AltermagnetSpinValve::crsb_cu_crsb();
//! let n = valve.harmonic_order();
//!
//! // Aligned layers: minimal magnetoresistance.
//! assert!(valve.gmr_ratio(0.0).abs() < 1e-9);
//!
//! // Channel-swap angle: maximal magnetoresistance.
//! let theta_max = PI / f64::from(n);
//! assert!((valve.gmr_ratio(theta_max) - valve.gmr_max).abs() < 1e-9);
//!
//! // Zero net magnetization at every configuration -- the whole point.
//! assert_eq!(valve.net_magnetization(), 0.0);
//! ```

use super::materials::Altermagnet;
use crate::error::{Error, Result};
use crate::material::multilayer::SpacerLayer;

/// Altermagnet-based spin valve exhibiting giant magnetoresistance (GMR)
/// without any ferromagnetic layer.
///
/// Two altermagnetic layers (each with zero net magnetization) sandwich a
/// non-magnetic metallic spacer. The device resistance depends on the
/// relative angle `theta` between the layers' Neel/crystal axes through the
/// shared harmonic order of their `AltermagneticSymmetry`, exactly as in a
/// conventional GMR spin valve -- but with no net moment anywhere in the
/// stack.
// Note: no serde support yet -- `SpacerLayer` (reused from
// `crate::material::multilayer`) does not implement Serialize/Deserialize,
// matching the same limitation on `MagneticMultilayer` itself.
#[derive(Debug, Clone)]
pub struct AltermagnetSpinValve {
    /// Bottom (reference) altermagnetic layer.
    pub bottom_layer: Altermagnet,
    /// Bottom layer thickness \[nm\].
    pub bottom_thickness: f64,
    /// Non-magnetic metallic spacer layer.
    pub spacer: SpacerLayer,
    /// Top (free) altermagnetic layer.
    pub top_layer: Altermagnet,
    /// Top layer thickness \[nm\].
    pub top_thickness: f64,
    /// Maximum GMR ratio (dimensionless fraction), reached at the
    /// channel-swap angle `theta = pi / n`.
    pub gmr_max: f64,
}

impl AltermagnetSpinValve {
    /// Create a new altermagnet spin valve.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - either layer thickness is not positive
    /// - `gmr_max` is not finite and positive
    /// - the two layers do not share the same harmonic order (`n` in
    ///   `AltermagneticSymmetry::harmonic_order`) -- the channel-swap angle
    ///   `pi/n` is only well-defined when both layers oscillate with the
    ///   same harmonic order
    pub fn new(
        bottom_layer: Altermagnet,
        bottom_thickness: f64,
        spacer: SpacerLayer,
        top_layer: Altermagnet,
        top_thickness: f64,
        gmr_max: f64,
    ) -> Result<Self> {
        if bottom_thickness <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "bottom_thickness".to_string(),
                reason: "Layer thickness must be positive".to_string(),
            });
        }
        if top_thickness <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "top_thickness".to_string(),
                reason: "Layer thickness must be positive".to_string(),
            });
        }
        if !gmr_max.is_finite() || gmr_max <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "gmr_max".to_string(),
                reason: "Maximum GMR ratio must be finite and positive".to_string(),
            });
        }
        let n_bottom = bottom_layer.symmetry.harmonic_order();
        let n_top = top_layer.symmetry.harmonic_order();
        if n_bottom != n_top {
            return Err(Error::InvalidParameter {
                param: "top_layer".to_string(),
                reason: format!(
                    "top layer harmonic order {n_top} must match bottom layer harmonic \
                     order {n_bottom} for a well-defined channel-swap angle"
                ),
            });
        }
        Ok(Self {
            bottom_layer,
            bottom_thickness,
            spacer,
            top_layer,
            top_thickness,
            gmr_max,
        })
    }

    /// Create a CrSb / Cu / CrSb altermagnet spin valve (g-wave, n = 4).
    ///
    /// CrSb is a well-characterized g-wave altermagnet with a high Neel
    /// temperature (~700 K), making it an attractive room-temperature
    /// candidate for GMR-without-ferromagnetism devices.
    pub fn crsb_cu_crsb() -> Self {
        Self {
            bottom_layer: Altermagnet::crsb(),
            bottom_thickness: 10.0,
            spacer: SpacerLayer::copper(2.0),
            top_layer: Altermagnet::crsb(),
            top_thickness: 5.0,
            gmr_max: 0.12,
        }
    }

    /// Create an MnTe-based altermagnet spin valve (d-wave, n = 2).
    ///
    /// MnTe was the first material in which altermagnetic spin splitting
    /// was directly observed via ARPES.
    pub fn mnte_based() -> Self {
        Self {
            bottom_layer: Altermagnet::mnte(),
            bottom_thickness: 8.0,
            spacer: SpacerLayer::copper(2.5),
            top_layer: Altermagnet::mnte(),
            top_thickness: 4.0,
            gmr_max: 0.08,
        }
    }

    /// Shared harmonic order `n` of the two altermagnetic layers.
    ///
    /// `n = 2` (d-wave), `4` (g-wave), or `6` (i-wave). Determines the
    /// channel-swap angle `pi / n` at which the GMR ratio is maximal.
    pub fn harmonic_order(&self) -> u32 {
        self.bottom_layer.symmetry.harmonic_order()
    }

    /// Compute the GMR ratio at relative Neel/crystal angle `theta` \[rad\].
    ///
    /// `GMR(theta) = gmr_max * (1 - cos(n*theta)) / 2`
    ///
    /// where `n` is the shared harmonic order. This is minimal (zero) at
    /// `theta = 0` (aligned layers) and maximal (`gmr_max`) at the
    /// channel-swap angle `theta = pi/n`.
    pub fn gmr_ratio(&self, theta: f64) -> f64 {
        let n = f64::from(self.harmonic_order());
        self.gmr_max * (1.0 - (n * theta).cos()) / 2.0
    }

    /// Compute the two-point resistance at relative angle `theta` \[rad\]
    /// over junction area `area` \[m²\].
    ///
    /// `R(theta, area) = R0 * (1 + GMR(theta)) / area`
    ///
    /// where `R0` is a fixed resistance-area product \[Ω·m²\], following the
    /// same convention as `MagneticMultilayer::resistance` in
    /// `src/material/multilayer.rs`. The resistance-area product
    /// `R(theta, area) * area` therefore depends only on `theta`, not on
    /// `area`.
    ///
    /// # Errors
    ///
    /// Returns an error if `area` is not finite and positive.
    pub fn resistance(&self, theta: f64, area: f64) -> Result<f64> {
        if !area.is_finite() || area <= 0.0 {
            return Err(Error::InvalidParameter {
                param: "area".to_string(),
                reason: "Junction area must be finite and positive".to_string(),
            });
        }
        let r0 = 1.0e-12; // Base resistance-area product [Ω·m²], matches MagneticMultilayer
        let mr = self.gmr_ratio(theta);
        Ok(r0 * (1.0 + mr) / area)
    }

    /// Compute the two-point conductance at relative angle `theta` \[rad\]
    /// over junction area `area` \[m²\].
    ///
    /// `G(theta, area) = 1 / R(theta, area)` \[S\] (Ohm's law), so
    /// `R(theta, area) * G(theta, area) == 1` exactly for any valid inputs.
    ///
    /// # Errors
    ///
    /// Returns an error if `area` is not finite and positive.
    pub fn conductance(&self, theta: f64, area: f64) -> Result<f64> {
        let r = self.resistance(theta, area)?;
        Ok(1.0 / r)
    }

    /// Net magnetization of the device \[A/m\].
    ///
    /// Always exactly zero: both layers are antiferromagnetically-ordered
    /// altermagnets whose sublattices cancel by construction, regardless of
    /// the relative angle `theta`, layer thicknesses, or spacer choice. This
    /// is the defining property that distinguishes this GMR-analog device
    /// from a conventional ferromagnetic spin valve.
    pub fn net_magnetization(&self) -> f64 {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::altermagnet::materials::AltermagneticSymmetry;
    use std::f64::consts::PI;

    #[test]
    fn test_crsb_cu_crsb_preset() {
        let valve = AltermagnetSpinValve::crsb_cu_crsb();
        assert_eq!(valve.bottom_layer.name, "CrSb");
        assert_eq!(valve.top_layer.name, "CrSb");
        assert_eq!(valve.harmonic_order(), 4);
        assert!(!valve.spacer.is_insulator);
        assert_eq!(valve.spacer.material, "Cu");
    }

    #[test]
    fn test_mnte_based_preset() {
        let valve = AltermagnetSpinValve::mnte_based();
        assert_eq!(valve.bottom_layer.name, "MnTe");
        assert_eq!(valve.top_layer.name, "MnTe");
        assert_eq!(valve.harmonic_order(), 2);
        assert!(!valve.spacer.is_insulator);
    }

    #[test]
    fn test_gmr_ratio_zero_when_aligned() {
        for valve in [
            AltermagnetSpinValve::crsb_cu_crsb(),
            AltermagnetSpinValve::mnte_based(),
        ] {
            let ratio = valve.gmr_ratio(0.0);
            assert!(
                ratio.abs() < 1e-9,
                "aligned layers should give ~zero GMR, got {ratio}"
            );
        }
    }

    #[test]
    fn test_gmr_ratio_maximal_at_channel_swap_angle() {
        for valve in [
            AltermagnetSpinValve::crsb_cu_crsb(),
            AltermagnetSpinValve::mnte_based(),
        ] {
            let n = f64::from(valve.harmonic_order());
            let theta_swap = PI / n;
            let ratio = valve.gmr_ratio(theta_swap);
            assert!(
                (ratio - valve.gmr_max).abs() < 1e-9,
                "GMR should reach gmr_max at the channel-swap angle, got {ratio} vs {}",
                valve.gmr_max
            );

            // It is truly a maximum: nearby angles give a smaller ratio.
            let just_below = valve.gmr_ratio(theta_swap - 1e-3);
            let just_above = valve.gmr_ratio(theta_swap + 1e-3);
            assert!(just_below < ratio + 1e-9);
            assert!(just_above < ratio + 1e-9);
        }
    }

    #[test]
    fn test_gmr_ratio_monotonic_on_channel_swap_interval() {
        let valve = AltermagnetSpinValve::crsb_cu_crsb();
        let n = f64::from(valve.harmonic_order());
        let theta_swap = PI / n;

        let mut prev = valve.gmr_ratio(0.0);
        for step in 1..=200 {
            let theta = theta_swap * f64::from(step) / 200.0;
            let ratio = valve.gmr_ratio(theta);
            assert!(
                ratio >= prev - 1e-12,
                "GMR ratio should increase monotonically on [0, pi/n]: \
                 theta={theta}, ratio={ratio}, prev={prev}"
            );
            prev = ratio;
        }
    }

    #[test]
    fn test_net_magnetization_always_zero() {
        let crsb_valve = AltermagnetSpinValve::crsb_cu_crsb();
        let mnte_valve = AltermagnetSpinValve::mnte_based();
        let custom_valve = AltermagnetSpinValve::new(
            Altermagnet::ruo2(),
            3.0,
            SpacerLayer::copper(1.5),
            Altermagnet::ruo2(),
            7.0,
            0.05,
        )
        .expect("valid custom spin valve");

        for valve in [crsb_valve, mnte_valve, custom_valve] {
            assert_eq!(
                valve.net_magnetization(),
                0.0,
                "GMR-without-ferromagnetism device must have exactly zero net moment"
            );
        }
    }

    #[test]
    fn test_resistance_conductance_consistency() {
        let valve = AltermagnetSpinValve::crsb_cu_crsb();
        let area = 1.0e-14;
        for &theta in &[0.0, 0.2, PI / 4.0, PI / 2.0, PI] {
            let r = valve.resistance(theta, area).expect("valid resistance");
            let g = valve.conductance(theta, area).expect("valid conductance");
            assert!(
                (r * g - 1.0).abs() < 1e-9,
                "R*G should equal 1 (Ohm's law): theta={theta}, R={r}, G={g}"
            );
        }

        // The resistance-area product depends only on theta, not on area,
        // matching the MagneticMultilayer convention (R = R0*(1+mr)/area).
        let n = f64::from(valve.harmonic_order());
        let theta = PI / n / 2.0;
        let ra_1 = valve.resistance(theta, 1.0e-14).expect("valid") * 1.0e-14;
        let ra_2 = valve.resistance(theta, 5.0e-13).expect("valid") * 5.0e-13;
        assert!(
            (ra_1 - ra_2).abs() / ra_1 < 1e-9,
            "resistance-area product should be independent of area: {ra_1} vs {ra_2}"
        );
    }

    #[test]
    fn test_resistance_channel_swap_greater_than_aligned() {
        let valve = AltermagnetSpinValve::crsb_cu_crsb();
        let area = 1.0e-14;
        let n = f64::from(valve.harmonic_order());

        let r_aligned = valve.resistance(0.0, area).expect("valid");
        let r_swapped = valve.resistance(PI / n, area).expect("valid");

        assert!(
            r_swapped > r_aligned,
            "channel-swapped configuration should have higher resistance"
        );
    }

    #[test]
    fn test_mismatched_harmonic_order_rejected() {
        let d_wave = Altermagnet::custom(
            "TestD",
            500.0,
            0.5,
            AltermagneticSymmetry::DWave,
            "test",
            1.0e5,
            4.0e-10,
        )
        .expect("valid d-wave material");
        let g_wave = Altermagnet::custom(
            "TestG",
            500.0,
            0.5,
            AltermagneticSymmetry::GWave,
            "test",
            1.0e5,
            4.0e-10,
        )
        .expect("valid g-wave material");

        let result =
            AltermagnetSpinValve::new(d_wave, 5.0, SpacerLayer::copper(2.0), g_wave, 5.0, 0.1);
        assert!(
            result.is_err(),
            "mismatched harmonic order must be rejected"
        );
    }

    #[test]
    fn test_invalid_thickness_rejected() {
        let result = AltermagnetSpinValve::new(
            Altermagnet::crsb(),
            0.0,
            SpacerLayer::copper(2.0),
            Altermagnet::crsb(),
            5.0,
            0.1,
        );
        assert!(result.is_err(), "zero bottom thickness must be rejected");

        let result = AltermagnetSpinValve::new(
            Altermagnet::crsb(),
            10.0,
            SpacerLayer::copper(2.0),
            Altermagnet::crsb(),
            -1.0,
            0.1,
        );
        assert!(result.is_err(), "negative top thickness must be rejected");
    }

    #[test]
    fn test_invalid_gmr_max_rejected() {
        let result = AltermagnetSpinValve::new(
            Altermagnet::crsb(),
            10.0,
            SpacerLayer::copper(2.0),
            Altermagnet::crsb(),
            5.0,
            0.0,
        );
        assert!(result.is_err(), "zero gmr_max must be rejected");

        let result = AltermagnetSpinValve::new(
            Altermagnet::crsb(),
            10.0,
            SpacerLayer::copper(2.0),
            Altermagnet::crsb(),
            5.0,
            f64::NAN,
        );
        assert!(result.is_err(), "non-finite gmr_max must be rejected");
    }

    #[test]
    fn test_invalid_area_rejected() {
        let valve = AltermagnetSpinValve::crsb_cu_crsb();
        assert!(valve.resistance(0.0, 0.0).is_err());
        assert!(valve.resistance(0.0, -1.0e-14).is_err());
        assert!(valve.conductance(0.0, 0.0).is_err());
    }
}
