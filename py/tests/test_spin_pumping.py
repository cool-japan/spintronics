"""Tests for the SpinPumpingSimulation high-level workflow binding
(src/python/simulation.rs): FMR precession -> spin pumping current -> ISHE
voltage detection.
"""
import math
import statistics

import pytest

import spintronics as sp

REQUIRED_RUN_KEYS = {
    "times",
    "mx",
    "my",
    "mz",
    "spin_current",
    "voltage",
    "peak_voltage",
    "avg_voltage",
    "peak_spin_current",
}


class TestConstructionAndDefaults:
    def test_default_state(self):
        sim = sp.SpinPumpingSimulation()
        assert sim.get_magnetization().to_tuple() == pytest.approx((1.0, 0.0, 0.0))

    def test_repr(self):
        text = repr(sp.SpinPumpingSimulation())
        assert text.startswith("SpinPumpingSimulation(")
        assert "YIG" in text
        assert "Pt" in text


class TestSetters:
    def test_set_magnetization_normalizes(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_magnetization(2.0, 0.0, 0.0)
        assert sim.get_magnetization().to_tuple() == pytest.approx((1.0, 0.0, 0.0))

    def test_set_fmr_conditions_resets_magnetization_along_x(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_magnetization(0.0, 1.0, 0.0)
        sim.set_fmr_conditions(9.65e9, 0.15)
        assert sim.get_magnetization().to_tuple() == pytest.approx((1.0, 0.0, 0.0))

    def test_set_sample_length_and_field_do_not_raise(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_sample_length(1.0e-2)
        sim.set_field(0.0, 0.0, 0.2)
        sim.set_magnetization(1.0, 0.0, 0.0)
        result = sim.run(1e-10, 20)
        assert len(result["times"]) == 21


class TestRun:
    def test_run_returns_all_expected_keys(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_fmr_conditions(9.65e9, 0.1)
        result = sim.run(1e-9, 200)
        assert set(result.keys()) == REQUIRED_RUN_KEYS

    def test_run_array_lengths_are_consistent(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_fmr_conditions(9.65e9, 0.1)
        n_steps = 200
        result = sim.run(1e-9, n_steps)
        expected_len = n_steps + 1
        assert len(result["times"]) == expected_len
        assert len(result["mx"]) == expected_len
        assert len(result["my"]) == expected_len
        assert len(result["mz"]) == expected_len
        assert len(result["spin_current"]) == expected_len
        assert len(result["voltage"]) == expected_len

    def test_run_starts_at_initial_conditions(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_fmr_conditions(9.65e9, 0.1)
        result = sim.run(1e-9, 100)
        assert result["times"][0] == 0.0
        assert (result["mx"][0], result["my"][0], result["mz"][0]) == pytest.approx(
            (1.0, 0.0, 0.0)
        )

    def test_run_final_time_matches_duration(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_fmr_conditions(9.65e9, 0.1)
        duration = 1e-9
        result = sim.run(duration, 200)
        assert result["times"][-1] == pytest.approx(duration)

    def test_magnetization_norm_conserved_throughout_run(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_fmr_conditions(9.65e9, 0.1)
        result = sim.run(1e-9, 300)
        for mx, my, mz in zip(result["mx"], result["my"], result["mz"]):
            norm = math.sqrt(mx * mx + my * my + mz * mz)
            assert norm == pytest.approx(1.0, abs=1e-9)

    def test_spin_current_and_voltage_are_nonnegative(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_fmr_conditions(9.65e9, 0.1)
        result = sim.run(1e-9, 200)
        assert all(js >= 0.0 for js in result["spin_current"])
        assert all(v >= 0.0 for v in result["voltage"])

    def test_peak_voltage_matches_max_of_series(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_fmr_conditions(9.65e9, 0.1)
        result = sim.run(1e-9, 200)
        assert result["peak_voltage"] == pytest.approx(max(result["voltage"]))

    def test_peak_spin_current_matches_max_of_series(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_fmr_conditions(9.65e9, 0.1)
        result = sim.run(1e-9, 200)
        assert result["peak_spin_current"] == pytest.approx(
            max(result["spin_current"])
        )

    def test_avg_voltage_matches_mean_of_series(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_fmr_conditions(9.65e9, 0.1)
        result = sim.run(1e-9, 200)
        assert result["avg_voltage"] == pytest.approx(
            statistics.mean(result["voltage"])
        )

    def test_larger_sample_length_gives_proportionally_larger_voltage(self):
        sim1 = sp.SpinPumpingSimulation()
        sim1.set_fmr_conditions(9.65e9, 0.1)
        sim1.set_sample_length(5e-3)
        result1 = sim1.run(1e-9, 200)

        sim2 = sp.SpinPumpingSimulation()
        sim2.set_fmr_conditions(9.65e9, 0.1)
        sim2.set_sample_length(1.0e-2)  # exactly double sim1's length
        result2 = sim2.run(1e-9, 200)

        assert result2["peak_voltage"] == pytest.approx(2.0 * result1["peak_voltage"])
        assert result2["avg_voltage"] == pytest.approx(2.0 * result1["avg_voltage"])

    def test_all_series_values_are_finite(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_fmr_conditions(9.65e9, 0.1)
        result = sim.run(1e-9, 200)
        for key in ("voltage", "spin_current", "mx", "my", "mz", "times"):
            assert all(math.isfinite(v) for v in result[key])

    def test_get_magnetization_before_run_matches_manual_set(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_magnetization(0.0, 1.0, 0.0)
        assert sim.get_magnetization().to_tuple() == pytest.approx((0.0, 1.0, 0.0))


class TestGetMagnetizationAfterRun:
    """Regression tests for ``get_magnetization()`` state after ``run()``.

    ``SpinPumpingSimulation.run()`` (src/python/simulation.rs) evolves a
    *local* copy of the magnetization (`let mut m = self.magnetization;`)
    internally and now writes the final value back to
    ``self.magnetization`` before returning the result dict, so
    ``get_magnetization()`` called after ``run()`` reports the trajectory's
    final state, matching the returned dict's mx/my/mz arrays (see TestRun
    tests above).

    This mirrors ``LlgSimulator``, whose ``evolve()`` also updates the
    internally tracked magnetization -- see
    test_llg_simulator.py::TestEvolve::
    test_evolve_updates_get_magnetization_to_final_state.
    """

    def test_get_magnetization_reflects_run_result(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_fmr_conditions(9.65e9, 0.1)
        result = sim.run(1e-9, 200)
        expected = (result["mx"][-1], result["my"][-1], result["mz"][-1])
        assert sim.get_magnetization().to_tuple() == pytest.approx(expected)

    def test_get_magnetization_changes_after_run(self):
        sim = sp.SpinPumpingSimulation()
        sim.set_fmr_conditions(9.65e9, 0.1)
        before = sim.get_magnetization().to_tuple()
        sim.run(1e-9, 200)
        after = sim.get_magnetization().to_tuple()
        distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(after, before)))
        assert distance > 1e-6
