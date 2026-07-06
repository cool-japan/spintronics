"""Bonus tests for the spin caloritronics bindings: OnsagerMatrix,
SpinCaloritronicsMaterial (src/python/caloritronics.rs).

See the module docstring in test_llb.py for a note on why these classes
are tested here despite not yet appearing in `__all__` / the `.pyi`
stub. Formulas and preset values were taken from
src/python/caloritronics.rs and cross-checked interactively.
"""
import pytest

import spintronics as sp


class TestOnsagerMatrixPresets:
    def test_yig_pt(self):
        mat = sp.OnsagerMatrix.yig_pt(300.0)
        assert mat.temperature == pytest.approx(300.0)
        assert mat.conductivity == pytest.approx(5.0e6)
        assert mat.seebeck == pytest.approx(-5.0e-6)
        assert mat.spin_seebeck == pytest.approx(1.0e-3)
        assert mat.hall_angle == pytest.approx(0.1)
        assert mat.thermal_conductivity == pytest.approx(46.0)

    def test_fe_pt_and_cofeb_pt_have_distinct_conductivity(self):
        fe_pt = sp.OnsagerMatrix.fe_pt(300.0)
        cofeb_pt = sp.OnsagerMatrix.cofeb_pt(300.0)
        assert fe_pt.conductivity == pytest.approx(1.0e7)
        assert cofeb_pt.conductivity == pytest.approx(8.0e6)

    def test_repr(self):
        assert repr(sp.OnsagerMatrix.yig_pt(300.0)).startswith("OnsagerMatrix(")


class TestOnsagerMatrixCustomConstruction:
    def test_custom_parameters(self):
        mat = sp.OnsagerMatrix(
            temperature=250.0,
            conductivity=1e6,
            seebeck=-5e-6,
            spin_seebeck=1e-3,
            hall_angle=0.01,
            thermal_conductivity=50.0,
        )
        assert mat.temperature == pytest.approx(250.0)
        assert mat.conductivity == pytest.approx(1e6)
        assert mat.seebeck == pytest.approx(-5e-6)
        assert mat.spin_seebeck == pytest.approx(1e-3)
        assert mat.hall_angle == pytest.approx(0.01)
        assert mat.thermal_conductivity == pytest.approx(50.0)


class TestOnsagerReciprocity:
    def test_analytically_constructed_matrices_satisfy_reciprocity(self):
        for mat in (
            sp.OnsagerMatrix.yig_pt(300.0),
            sp.OnsagerMatrix.fe_pt(300.0),
            sp.OnsagerMatrix.cofeb_pt(300.0),
        ):
            assert mat.reciprocity_error() < 1e-10


class TestOnsagerTransportFormulas:
    def test_spin_current_from_grad_t(self):
        mat = sp.OnsagerMatrix.yig_pt(300.0)
        js = mat.spin_current_from_grad_t((1000.0, 0.0, 0.0))
        # j_s = S_s * sigma * grad_T = 1e-3 * 5e6 * 1000 = 5e6
        assert js == pytest.approx((5.0e6, 0.0, 0.0))

    def test_heat_current_from_spin_current(self):
        mat = sp.OnsagerMatrix.yig_pt(300.0)
        jq = mat.heat_current_from_spin_current((5.0e6, 0.0, 0.0))
        # j_Q = T * S_s * j_s = 300 * 1e-3 * 5e6 = 1.5e6
        assert jq == pytest.approx((1.5e6, 0.0, 0.0))

    def test_nernst_voltage(self):
        mat = sp.OnsagerMatrix.yig_pt(300.0)
        # nu = -S_e * theta_H * |grad_T| = 5e-6 * 0.1 * 1000 = 5e-4
        assert mat.nernst_voltage(1000.0) == pytest.approx(5.0e-4)

    def test_all_currents_returns_expected_keys_and_values(self):
        mat = sp.OnsagerMatrix.yig_pt(300.0)
        currents = mat.all_currents((1000.0, 0.0, 0.0), (0.0, 0.0, 0.0))
        assert set(currents.keys()) == {
            "charge_current",
            "spin_current",
            "heat_current",
        }
        assert currents["charge_current"] == pytest.approx((-25000.0, 0.0, 0.0))
        assert currents["spin_current"] == pytest.approx((5.0e6, 0.0, 0.0))
        assert currents["heat_current"] == pytest.approx((-46000.0, 0.0, 0.0), abs=1e-6)


class TestSpinCaloritronicsMaterial:
    def test_presets_have_expected_names(self):
        assert sp.SpinCaloritronicsMaterial.yig_pt(300.0).name == "YIG/Pt"
        assert sp.SpinCaloritronicsMaterial.fe_pt(300.0).name == "Fe/Pt"
        assert sp.SpinCaloritronicsMaterial.cofeb_pt(300.0).name == "CoFeB/Pt"

    def test_repr(self):
        text = repr(sp.SpinCaloritronicsMaterial.yig_pt(300.0))
        assert text.startswith("SpinCaloritronicsMaterial(")

    def test_onsager_property_round_trips_temperature(self):
        scm = sp.SpinCaloritronicsMaterial.yig_pt(300.0)
        assert scm.onsager.temperature == pytest.approx(300.0)

    def test_from_onsager(self):
        onsager = sp.OnsagerMatrix.yig_pt(300.0)
        scm = sp.SpinCaloritronicsMaterial.from_onsager(onsager)
        assert scm.onsager.temperature == pytest.approx(300.0)

    def test_compute_all_returns_expected_keys(self):
        scm = sp.SpinCaloritronicsMaterial.yig_pt(300.0)
        result = scm.compute_all(grad_t=(1000.0, 0.0, 0.0), j_spin=(0.0, 0.0, 1e3))
        assert set(result.keys()) == {
            "spin_seebeck_current",
            "peltier_heat",
            "nernst_voltage",
            "spin_nernst_current",
            "reciprocity_satisfied",
            "total_heat_current",
        }

    def test_compute_all_values(self):
        scm = sp.SpinCaloritronicsMaterial.yig_pt(300.0)
        result = scm.compute_all(grad_t=(1000.0, 0.0, 0.0), j_spin=(0.0, 0.0, 1e3))
        assert result["spin_seebeck_current"] == pytest.approx((5.0e6, 0.0, 0.0))
        assert result["peltier_heat"] == pytest.approx(300.0)
        assert result["nernst_voltage"] == pytest.approx(5.0e-4)
        assert result["spin_nernst_current"] == pytest.approx((5.0e5, 0.0, 0.0))
        assert result["reciprocity_satisfied"] is True
        assert result["total_heat_current"] == pytest.approx(
            (-46000.0, 0.0, 300.0), abs=1e-6
        )
