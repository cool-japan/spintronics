"""Bonus tests for the Landau-Lifshitz-Bloch (LLB) finite-temperature
dynamics bindings: LlbMaterial, LlbSolver (src/python/llb.rs).

These classes are registered on the compiled extension module and are
reachable via `import spintronics; spintronics.LlbMaterial` (see
`from .spintronics import *` in python/spintronics/__init__.py and the
`m.add_class::<PyLlbMaterial>()` etc. registrations in src/python/mod.rs),
but are not yet listed in `__all__` / the `.pyi` type stub as of this
writing. Signatures and preset values were taken directly from the
authoritative Rust source (src/python/llb.rs, src/dynamics/llb.rs) and
cross-checked interactively before writing these assertions.
"""
import math

import pytest

import spintronics as sp


class TestLlbMaterialPresets:
    def test_iron(self):
        iron = sp.LlbMaterial.iron()
        assert iron.curie_temp == pytest.approx(1043.0)
        assert iron.alpha == pytest.approx(0.01)
        assert iron.spin_s == pytest.approx(1.0)
        assert iron.ms_0 == pytest.approx(1.71e6)

    def test_nickel(self):
        nickel = sp.LlbMaterial.nickel()
        assert nickel.curie_temp == pytest.approx(631.0)
        assert nickel.alpha == pytest.approx(0.064)
        assert nickel.spin_s == pytest.approx(0.3)
        assert nickel.ms_0 == pytest.approx(4.84e5)

    def test_cofeb(self):
        cofeb = sp.LlbMaterial.cofeb()
        assert cofeb.curie_temp == pytest.approx(1000.0)
        assert cofeb.alpha == pytest.approx(0.005)
        assert cofeb.spin_s == pytest.approx(0.5)
        assert cofeb.ms_0 == pytest.approx(1.2e6)

    def test_repr(self):
        assert repr(sp.LlbMaterial.iron()).startswith("LlbMaterial(")


class TestLlbMaterialCustomConstruction:
    def test_custom_parameters(self):
        mat = sp.LlbMaterial(curie_temp=800.0, alpha=0.02, spin_s=0.5, ms_0=1e6)
        assert mat.curie_temp == pytest.approx(800.0)
        assert mat.alpha == pytest.approx(0.02)
        assert mat.spin_s == pytest.approx(0.5)
        assert mat.ms_0 == pytest.approx(1e6)


class TestEquilibriumMagnetization:
    def test_below_curie_temperature_is_between_zero_and_one(self):
        iron = sp.LlbMaterial.iron()
        m_e = iron.equilibrium_magnetization(300.0)
        assert 0.0 < m_e < 1.0

    def test_vanishes_at_and_above_curie_temperature(self):
        iron = sp.LlbMaterial.iron()
        assert iron.equilibrium_magnetization(iron.curie_temp) == 0.0
        assert iron.equilibrium_magnetization(iron.curie_temp + 100.0) == 0.0

    def test_vanishes_at_and_below_absolute_zero(self):
        iron = sp.LlbMaterial.iron()
        assert iron.equilibrium_magnetization(0.0) == 0.0
        assert iron.equilibrium_magnetization(-10.0) == 0.0

    def test_decreases_with_increasing_temperature(self):
        iron = sp.LlbMaterial.iron()
        m_low = iron.equilibrium_magnetization(100.0)
        m_high = iron.equilibrium_magnetization(900.0)
        assert m_low > m_high


class TestDamping:
    def test_alpha_parallel_and_perp_are_positive_below_curie(self):
        iron = sp.LlbMaterial.iron()
        assert iron.alpha_parallel(300.0) > 0.0
        assert iron.alpha_perp(300.0) > 0.0


class TestLlbSolver:
    def test_construction(self):
        iron = sp.LlbMaterial.iron()
        solver = sp.LlbSolver(iron, dt=1e-14, temperature=300.0, h_ext=(0.0, 0.0, 1.0))
        assert solver.dt == pytest.approx(1e-14)
        assert solver.temperature == pytest.approx(300.0)
        assert solver.h_ext == pytest.approx((0.0, 0.0, 1.0))
        assert solver.gamma == pytest.approx(sp.GAMMA)

    def test_repr(self):
        iron = sp.LlbMaterial.iron()
        solver = sp.LlbSolver(iron, dt=1e-14, temperature=300.0, h_ext=(0.0, 0.0, 1.0))
        assert repr(solver).startswith("LlbSolver(")

    def test_step_returns_three_finite_components(self):
        iron = sp.LlbMaterial.iron()
        solver = sp.LlbSolver(iron, dt=1e-14, temperature=300.0, h_ext=(0.0, 0.0, 1.0))
        m_new = solver.step((0.8, 0.1, 0.0))
        assert len(m_new) == 3
        assert all(math.isfinite(c) for c in m_new)

    def test_run_returns_expected_keys_and_lengths(self):
        iron = sp.LlbMaterial.iron()
        solver = sp.LlbSolver(iron, dt=1e-14, temperature=300.0, h_ext=(0.0, 0.0, 1.0))
        result = solver.run(m0=(0.8, 0.1, 0.0), num_steps=1000, record_every=10)
        assert set(result.keys()) == {
            "mx",
            "my",
            "mz",
            "m_magnitude",
            "time",
            "equilibrium_m",
        }
        # 1000 steps recorded every 10th step, plus the initial snapshot.
        assert len(result["mx"]) == 101
        assert len(result["time"]) == 101
        assert result["time"][0] == pytest.approx(0.0)
        assert result["time"][-1] == pytest.approx(1000 * 1e-14)
        assert result["equilibrium_m"] == pytest.approx(
            iron.equilibrium_magnetization(300.0)
        )

    def test_run_initial_magnitude_matches_m0(self):
        iron = sp.LlbMaterial.iron()
        solver = sp.LlbSolver(iron, dt=1e-14, temperature=300.0, h_ext=(0.0, 0.0, 1.0))
        m0 = (0.8, 0.1, 0.0)
        result = solver.run(m0=m0, num_steps=100, record_every=10)
        expected_mag = math.sqrt(m0[0] ** 2 + m0[1] ** 2 + m0[2] ** 2)
        assert result["m_magnitude"][0] == pytest.approx(expected_mag)
