"""Tests for the InverseSpinHall (ISHE) binding (src/python/effects.rs),
cross-checked against the formulas in src/effect/ishe.rs:

    E = rho * theta_sh * (js_flow x js_polarization)      [V/m]
    voltage = |E| * strip_width                            [V]
    efficiency = rho * |theta_sh|
"""
import pytest

import spintronics as sp


class TestPresets:
    def test_platinum(self):
        pt = sp.InverseSpinHall.platinum()
        assert pt.theta_sh == pytest.approx(0.08)
        assert pt.rho == pytest.approx(2.0e-7)

    def test_tantalum(self):
        ta = sp.InverseSpinHall.tantalum()
        assert ta.theta_sh == pytest.approx(0.12)
        assert ta.rho == pytest.approx(1.8e-7)

    def test_tungsten(self):
        w = sp.InverseSpinHall.tungsten()
        assert w.theta_sh == pytest.approx(-0.30)
        assert w.rho == pytest.approx(5.0e-8)

    def test_tungsten_has_negative_spin_hall_angle(self):
        # Tungsten's giant negative spin Hall angle has opposite sign to
        # both Pt and Ta.
        assert sp.InverseSpinHall.tungsten().theta_sh < 0.0
        assert sp.InverseSpinHall.platinum().theta_sh > 0.0
        assert sp.InverseSpinHall.tantalum().theta_sh > 0.0


class TestCustomConstruction:
    def test_custom_parameters(self):
        ishe = sp.InverseSpinHall(theta_sh=0.05, rho=1.0e-7)
        assert ishe.theta_sh == pytest.approx(0.05)
        assert ishe.rho == pytest.approx(1.0e-7)

    def test_repr(self, platinum):
        assert repr(platinum).startswith("InverseSpinHall(")


class TestConvert:
    def test_perpendicular_flow_and_polarization(self, platinum):
        js_flow = sp.Vector3(0.0, 0.0, 1.0)
        js_pol = sp.Vector3(1.0e8, 0.0, 0.0)
        e_field = platinum.convert(js_flow, js_pol)
        # E = rho * theta_sh * (js_flow x js_pol)
        #   = 2e-7 * 0.08 * ((0,0,1) x (1e8,0,0))
        #   = 1.6e-8 * (0, 1e8, 0) = (0, 1.6, 0)
        assert e_field.to_tuple() == pytest.approx((0.0, 1.6, 0.0))

    def test_field_is_perpendicular_to_both_inputs(self, platinum):
        js_flow = sp.Vector3(1.0, 0.0, 0.0)
        js_pol = sp.Vector3(0.0, 1.0, 0.0)
        e_field = platinum.convert(js_flow, js_pol)
        assert e_field.dot(js_flow) == pytest.approx(0.0, abs=1e-12)
        assert e_field.dot(js_pol) == pytest.approx(0.0, abs=1e-12)

    def test_zero_current_gives_zero_field(self, platinum):
        js_flow = sp.Vector3.zero()
        js_pol = sp.Vector3(0.0, 1.0, 0.0)
        e_field = platinum.convert(js_flow, js_pol)
        assert e_field.magnitude() == pytest.approx(0.0, abs=1e-30)

    def test_parallel_vectors_give_zero_field(self, platinum):
        v = sp.Vector3(1.0, 0.0, 0.0)
        e_field = platinum.convert(v, v)
        assert e_field.magnitude() == pytest.approx(0.0, abs=1e-30)

    def test_tungsten_field_has_opposite_sign_to_platinum(self):
        js_flow = sp.Vector3(0.0, 0.0, 1.0)
        js_pol = sp.Vector3(1.0e8, 0.0, 0.0)
        pt_field = sp.InverseSpinHall.platinum().convert(js_flow, js_pol)
        w_field = sp.InverseSpinHall.tungsten().convert(js_flow, js_pol)
        assert pt_field.y > 0.0
        assert w_field.y < 0.0


class TestVoltage:
    def test_voltage_matches_field_magnitude_times_width(self, platinum):
        js_flow = sp.Vector3(0.0, 0.0, 1.0)
        js_pol = sp.Vector3(1.0e8, 0.0, 0.0)
        width = 1.0e-4
        voltage = platinum.voltage(js_flow, js_pol, width)
        e_field = platinum.convert(js_flow, js_pol)
        assert voltage == pytest.approx(e_field.magnitude() * width)

    def test_voltage_is_nonnegative(self):
        w = sp.InverseSpinHall.tungsten()
        js_flow = sp.Vector3(0.0, 0.0, 1.0)
        js_pol = sp.Vector3(1.0e8, 0.0, 0.0)
        voltage = w.voltage(js_flow, js_pol, 1e-4)
        assert voltage >= 0.0

    def test_voltage_scales_linearly_with_width(self, platinum):
        js_flow = sp.Vector3(0.0, 0.0, 1.0)
        js_pol = sp.Vector3(1.0e8, 0.0, 0.0)
        v1 = platinum.voltage(js_flow, js_pol, 1.0e-4)
        v2 = platinum.voltage(js_flow, js_pol, 2.0e-4)
        assert v2 == pytest.approx(2.0 * v1)


class TestEfficiency:
    def test_efficiency_formula(self, platinum):
        assert platinum.efficiency() == pytest.approx(platinum.rho * abs(platinum.theta_sh))

    def test_efficiency_is_nonnegative_even_for_negative_theta(self):
        w = sp.InverseSpinHall.tungsten()
        assert w.theta_sh < 0.0
        assert w.efficiency() > 0.0
        assert w.efficiency() == pytest.approx(w.rho * abs(w.theta_sh))
