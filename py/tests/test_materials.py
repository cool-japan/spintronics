"""Tests for the Ferromagnet and SpinInterface material bindings
(src/python/materials.rs), cross-checked against the underlying Rust
presets in src/material/ferromagnet.rs and src/material/interface.rs.
"""
import pytest

import spintronics as sp

# ---------------------------------------------------------------------------
# Ferromagnet
# ---------------------------------------------------------------------------


class TestFerromagnetPresets:
    def test_yig(self):
        yig = sp.Ferromagnet.yig()
        assert yig.alpha == pytest.approx(0.0001)
        assert yig.ms == pytest.approx(1.4e5)
        assert yig.anisotropy_k == pytest.approx(0.0)
        assert yig.easy_axis.to_tuple() == pytest.approx((0.0, 0.0, 1.0))
        assert yig.exchange_a == pytest.approx(3.7e-12)

    def test_yig_has_ultra_low_damping(self):
        # YIG's extremely low Gilbert damping is what makes it ideal for
        # spin pumping / magnon transport experiments.
        assert sp.Ferromagnet.yig().alpha < sp.Ferromagnet.permalloy().alpha

    def test_permalloy(self):
        py = sp.Ferromagnet.permalloy()
        assert py.alpha == pytest.approx(0.01)
        assert py.ms == pytest.approx(8.0e5)
        assert py.anisotropy_k == pytest.approx(0.0)
        assert py.easy_axis.to_tuple() == pytest.approx((1.0, 0.0, 0.0))
        assert py.exchange_a == pytest.approx(1.3e-11)

    def test_cofe(self):
        cofe = sp.Ferromagnet.cofe()
        assert cofe.alpha == pytest.approx(0.003)
        assert cofe.ms == pytest.approx(1.95e6)
        assert cofe.anisotropy_k == pytest.approx(1.0e4)
        assert cofe.easy_axis.to_tuple() == pytest.approx((0.0, 0.0, 1.0))
        assert cofe.exchange_a == pytest.approx(2.5e-11)

    def test_cofeb(self):
        cofeb = sp.Ferromagnet.cofeb()
        assert cofeb.alpha == pytest.approx(0.004)
        assert cofeb.ms == pytest.approx(1.0e6)
        assert cofeb.anisotropy_k == pytest.approx(1.0e5)
        assert cofeb.easy_axis.to_tuple() == pytest.approx((0.0, 0.0, 1.0))
        assert cofeb.exchange_a == pytest.approx(2.0e-11)

    def test_iron(self):
        iron = sp.Ferromagnet.iron()
        assert iron.alpha == pytest.approx(0.002)
        assert iron.ms == pytest.approx(1.71e6)
        assert iron.anisotropy_k == pytest.approx(4.8e4)
        assert iron.easy_axis.to_tuple() == pytest.approx((1.0, 0.0, 0.0))
        assert iron.exchange_a == pytest.approx(2.1e-11)

    def test_cobalt(self):
        cobalt = sp.Ferromagnet.cobalt()
        assert cobalt.alpha == pytest.approx(0.005)
        assert cobalt.ms == pytest.approx(1.4e6)
        assert cobalt.anisotropy_k == pytest.approx(5.0e5)
        assert cobalt.easy_axis.to_tuple() == pytest.approx((0.0, 0.0, 1.0))
        assert cobalt.exchange_a == pytest.approx(3.0e-11)

    def test_nickel(self):
        nickel = sp.Ferromagnet.nickel()
        assert nickel.alpha == pytest.approx(0.045)
        assert nickel.ms == pytest.approx(4.8e5)
        assert nickel.anisotropy_k == pytest.approx(-4.5e3)
        assert nickel.exchange_a == pytest.approx(0.9e-11)

    def test_nickel_has_highest_damping_among_presets(self):
        presets = {
            "yig": sp.Ferromagnet.yig(),
            "permalloy": sp.Ferromagnet.permalloy(),
            "cofe": sp.Ferromagnet.cofe(),
            "cofeb": sp.Ferromagnet.cofeb(),
            "iron": sp.Ferromagnet.iron(),
            "cobalt": sp.Ferromagnet.cobalt(),
            "nickel": sp.Ferromagnet.nickel(),
        }
        highest = max(presets, key=lambda name: presets[name].alpha)
        assert highest == "nickel"


class TestFerromagnetCustomConstruction:
    def test_required_parameters_with_defaults(self):
        fm = sp.Ferromagnet(alpha=0.02, ms=5.0e5)
        assert fm.alpha == pytest.approx(0.02)
        assert fm.ms == pytest.approx(5.0e5)
        # Defaults per the #[pyo3(signature = ...)] binding.
        assert fm.anisotropy_k == 0.0
        assert fm.easy_axis.to_tuple() == pytest.approx((0.0, 0.0, 1.0))
        assert fm.exchange_a == pytest.approx(1e-11)

    def test_all_parameters_positional(self):
        fm = sp.Ferromagnet(0.02, 5.0e5, 1.0e4, (1.0, 0.0, 0.0), 2.0e-11)
        assert fm.alpha == pytest.approx(0.02)
        assert fm.ms == pytest.approx(5.0e5)
        assert fm.anisotropy_k == pytest.approx(1.0e4)
        assert fm.easy_axis.to_tuple() == pytest.approx((1.0, 0.0, 0.0))
        assert fm.exchange_a == pytest.approx(2.0e-11)

    def test_all_parameters_keyword(self):
        fm = sp.Ferromagnet(
            alpha=0.02,
            ms=5.0e5,
            anisotropy_k=1.0e4,
            easy_axis=(1.0, 0.0, 0.0),
            exchange_a=2.0e-11,
        )
        assert fm.anisotropy_k == pytest.approx(1.0e4)
        assert fm.exchange_a == pytest.approx(2.0e-11)

    def test_easy_axis_is_normalized(self):
        fm = sp.Ferromagnet(alpha=0.01, ms=1e5, easy_axis=(2.0, 0.0, 0.0))
        assert fm.easy_axis.to_tuple() == pytest.approx((1.0, 0.0, 0.0))
        assert fm.easy_axis.magnitude() == pytest.approx(1.0)

    def test_easy_axis_diagonal_normalization(self):
        fm = sp.Ferromagnet(alpha=0.01, ms=1e5, easy_axis=(1.0, 1.0, 0.0))
        expected = 1.0 / (2.0**0.5)
        assert fm.easy_axis.to_tuple() == pytest.approx((expected, expected, 0.0))

    def test_repr(self, yig):
        text = repr(yig)
        assert text.startswith("Ferromagnet(")
        assert "alpha=" in text
        assert "ms=" in text


# ---------------------------------------------------------------------------
# SpinInterface
# ---------------------------------------------------------------------------


class TestSpinInterfacePresets:
    def test_yig_pt(self):
        interface = sp.SpinInterface.yig_pt()
        assert interface.g_r == pytest.approx(1.0e19)
        assert interface.g_i == pytest.approx(0.0)
        assert interface.normal.to_tuple() == pytest.approx((0.0, 1.0, 0.0))
        assert interface.area == pytest.approx(1.0e-12)

    def test_py_pt(self):
        interface = sp.SpinInterface.py_pt()
        assert interface.g_r == pytest.approx(5.0e19)
        assert interface.g_i == pytest.approx(0.0)
        assert interface.normal.to_tuple() == pytest.approx((0.0, 1.0, 0.0))
        assert interface.area == pytest.approx(1.0e-12)

    def test_py_pt_has_higher_conductance_than_yig_pt(self):
        # Metallic Permalloy has a much higher spin-mixing conductance at
        # a Pt interface than the insulating YIG.
        assert sp.SpinInterface.py_pt().g_r > sp.SpinInterface.yig_pt().g_r


class TestSpinInterfaceCustomConstruction:
    def test_required_and_default_parameters(self):
        interface = sp.SpinInterface(g_r=2.0e19)
        assert interface.g_r == pytest.approx(2.0e19)
        assert interface.g_i == 0.0
        assert interface.normal.to_tuple() == pytest.approx((0.0, 1.0, 0.0))
        assert interface.area == pytest.approx(1e-12)

    def test_all_parameters(self):
        interface = sp.SpinInterface(
            g_r=2.0e19, g_i=1.0e18, normal=(1.0, 0.0, 0.0), area=5.0e-13
        )
        assert interface.g_r == pytest.approx(2.0e19)
        assert interface.g_i == pytest.approx(1.0e18)
        assert interface.normal.to_tuple() == pytest.approx((1.0, 0.0, 0.0))
        assert interface.area == pytest.approx(5.0e-13)

    def test_normal_is_not_normalized(self):
        # Unlike Ferromagnet.easy_axis, SpinInterface does NOT normalize the
        # `normal` argument (verified against src/python/materials.rs::
        # PySpinInterface::new, which stores Vector3::new(...) directly with
        # no .normalize() call). This is consistent with the .pyi docstring,
        # which -- unlike Ferromagnet.easy_axis's "normalised Vector3" -- makes
        # no normalization claim for SpinInterface.normal.
        interface = sp.SpinInterface(g_r=1e19, normal=(2.0, 0.0, 0.0))
        assert interface.normal.to_tuple() == pytest.approx((2.0, 0.0, 0.0))
        assert interface.normal.magnitude() == pytest.approx(2.0)

    def test_repr(self, yig_pt_interface):
        text = repr(yig_pt_interface)
        assert text.startswith("SpinInterface(")
        assert "g_r=" in text
