"""Tests for the module-level physical constants exposed by the spintronics
extension module (see src/python/mod.rs and src/constants.rs).
"""
import math

import pytest

import spintronics as sp


class TestFundamentalConstants:
    def test_hbar(self):
        # Reduced Planck constant [J.s] (CODATA 2018 exact value).
        assert sp.HBAR == pytest.approx(1.054571817e-34, rel=1e-12)

    def test_e_charge(self):
        # Elementary charge [C] (exact by SI definition since 2019).
        assert sp.E_CHARGE == pytest.approx(1.602176634e-19, rel=1e-12)

    def test_kb(self):
        # Boltzmann constant [J/K] (exact by SI definition since 2019).
        assert sp.KB == pytest.approx(1.380649e-23, rel=1e-12)

    def test_mu_b(self):
        # Bohr magneton [J/T].
        assert sp.MU_B == pytest.approx(9.2740100783e-24, rel=1e-9)

    def test_gamma(self):
        # Gyromagnetic ratio of the electron [rad/(s.T)].
        assert sp.GAMMA == pytest.approx(1.76085963023e11, rel=1e-9)

    @pytest.mark.parametrize("name", ["HBAR", "GAMMA", "E_CHARGE", "MU_B", "KB"])
    def test_constants_are_positive_floats(self, name):
        value = getattr(sp, name)
        assert isinstance(value, float)
        assert value > 0.0


class TestDerivedRelationships:
    """Sanity-check well-known physical relationships between the constants,
    rather than just re-asserting the literal values above."""

    def test_thermal_voltage_at_room_temperature(self):
        # V_T = k_B*T/e at T = 300 K is the textbook ~25.85 mV thermal voltage.
        v_t = sp.KB * 300.0 / sp.E_CHARGE
        assert v_t == pytest.approx(0.02585, rel=1e-3)

    def test_electron_gyromagnetic_ratio_over_two_pi(self):
        # gamma / 2*pi is the well-known ~28 GHz/T electron spin resonance
        # frequency-to-field conversion factor.
        freq_per_tesla_ghz = sp.GAMMA / (2.0 * math.pi) / 1e9
        assert freq_per_tesla_ghz == pytest.approx(28.0, rel=0.02)

    def test_gamma_mu_b_hbar_lande_relation(self):
        # gamma = g * mu_B / hbar, with the free-electron g-factor ~2.0023.
        g_factor = sp.GAMMA * sp.HBAR / sp.MU_B
        assert g_factor == pytest.approx(2.0023, rel=1e-3)
