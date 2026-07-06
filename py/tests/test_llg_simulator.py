"""Tests for the LlgSimulator (Landau-Lifshitz-Gilbert) binding
(src/python/dynamics.rs).
"""
import math

import pytest

import spintronics as sp


class TestConstructionAndDefaults:
    def test_default_state(self, yig):
        sim = sp.LlgSimulator(yig)
        assert sim.get_magnetization().to_tuple() == pytest.approx((1.0, 0.0, 0.0))
        assert sim.get_external_field().to_tuple() == pytest.approx((0.0, 0.0, 0.0))
        assert sim.time == 0.0

    def test_repr(self, yig):
        assert repr(sp.LlgSimulator(yig)).startswith("LlgSimulator(")


class TestSettersAndGetters:
    def test_set_and_get_magnetization(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_magnetization(0.0, 1.0, 0.0)
        assert sim.get_magnetization().to_tuple() == pytest.approx((0.0, 1.0, 0.0))

    def test_set_magnetization_normalizes(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_magnetization(2.0, 0.0, 0.0)
        m = sim.get_magnetization()
        assert m.to_tuple() == pytest.approx((1.0, 0.0, 0.0))
        assert m.magnitude() == pytest.approx(1.0)

    def test_set_and_get_external_field(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_external_field(0.0, 0.0, 0.2)
        assert sim.get_external_field().to_tuple() == pytest.approx((0.0, 0.0, 0.2))

    def test_reset_time(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_external_field(0.0, 0.0, 0.1)
        sim.step_rk4(1e-12)
        assert sim.time > 0.0
        sim.reset_time()
        assert sim.time == 0.0


class TestDmDt:
    def test_dm_dt_perpendicular_to_magnetization(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_magnetization(1.0, 0.0, 0.0)
        sim.set_external_field(0.0, 0.0, 0.1)
        dm_dt = sim.dm_dt()
        assert dm_dt.dot(sim.get_magnetization()) == pytest.approx(0.0, abs=1e-8)

    def test_zero_field_gives_zero_dm_dt(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_magnetization(1.0, 0.0, 0.0)
        sim.set_external_field(0.0, 0.0, 0.0)
        assert sim.dm_dt().to_tuple() == pytest.approx((0.0, 0.0, 0.0), abs=1e-15)


class TestPrecession:
    def test_precession_frequency_matches_larmor_formula(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_external_field(0.0, 0.0, 0.1)
        assert sim.precession_frequency() == pytest.approx(sp.GAMMA * 0.1)

    def test_precession_period_matches_frequency(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_external_field(0.0, 0.0, 0.1)
        expected_period = 2.0 * math.pi / (sp.GAMMA * 0.1)
        assert sim.precession_period() == pytest.approx(expected_period)

    def test_zero_field_gives_infinite_period(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_external_field(0.0, 0.0, 0.0)
        assert sim.precession_frequency() == 0.0
        assert math.isinf(sim.precession_period())


class TestSteppingIntegrators:
    def test_step_rk4_advances_time(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_external_field(0.0, 0.0, 0.1)
        dt = 1e-12
        sim.step_rk4(dt)
        assert sim.time == pytest.approx(dt)

    def test_step_euler_advances_time(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_external_field(0.0, 0.0, 0.1)
        dt = 1e-12
        sim.step_euler(dt)
        assert sim.time == pytest.approx(dt)

    def test_step_rk4_conserves_norm(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_magnetization(1.0, 0.0, 0.0)
        sim.set_external_field(0.0, 0.0, 0.1)
        for _ in range(200):
            sim.step_rk4(1e-12)
        assert sim.get_magnetization().magnitude() == pytest.approx(1.0, abs=1e-9)

    def test_step_euler_conserves_norm_via_renormalization(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_magnetization(1.0, 0.0, 0.0)
        sim.set_external_field(0.0, 0.0, 0.1)
        for _ in range(200):
            sim.step_euler(1e-13)
        assert sim.get_magnetization().magnitude() == pytest.approx(1.0, abs=1e-9)


class TestEvolve:
    def test_evolve_returns_n_steps_plus_one_points(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_magnetization(1.0, 0.0, 0.0)
        sim.set_external_field(0.0, 0.0, 0.1)
        traj = sim.evolve(1e-9, 100)
        assert len(traj) == 101

    def test_evolve_starts_at_initial_state(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_magnetization(1.0, 0.0, 0.0)
        sim.set_external_field(0.0, 0.0, 0.1)
        traj = sim.evolve(1e-9, 50)
        t0, mx0, my0, mz0 = traj[0]
        assert t0 == 0.0
        assert (mx0, my0, mz0) == pytest.approx((1.0, 0.0, 0.0))

    def test_evolve_default_method_is_rk4(self, yig):
        sim_default = sp.LlgSimulator(yig)
        sim_default.set_magnetization(1.0, 0.0, 0.0)
        sim_default.set_external_field(0.0, 0.0, 0.1)
        traj_default = sim_default.evolve(1e-9, 100)

        sim_rk4 = sp.LlgSimulator(yig)
        sim_rk4.set_magnetization(1.0, 0.0, 0.0)
        sim_rk4.set_external_field(0.0, 0.0, 0.1)
        traj_rk4 = sim_rk4.evolve(1e-9, 100, method="rk4")

        assert traj_default[-1] == pytest.approx(traj_rk4[-1])

    def test_evolve_advances_final_time_to_duration(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_external_field(0.0, 0.0, 0.1)
        duration = 1e-9
        traj = sim.evolve(duration, 100)
        assert traj[-1][0] == pytest.approx(duration)
        assert sim.time == pytest.approx(duration)

    def test_evolve_updates_get_magnetization_to_final_state(self, yig):
        # get_magnetization() should reflect the state after evolve(), the
        # same way it does after step_rk4()/step_euler(). Contrast this with
        # SpinPumpingSimulation.run(), which currently does NOT update its
        # tracked magnetization (see test_spin_pumping.py::
        # TestGetMagnetizationAfterRun).
        sim = sp.LlgSimulator(yig)
        sim.set_magnetization(1.0, 0.0, 0.0)
        sim.set_external_field(0.0, 0.0, 0.1)
        traj = sim.evolve(1e-9, 100)
        _, mx, my, mz = traj[-1]
        assert sim.get_magnetization().to_tuple() == pytest.approx((mx, my, mz))

    def test_rk4_and_euler_diverge_for_coarse_steps(self, yig):
        # rk4 (4th order) and euler (1st order) should give noticeably
        # different trajectories once dt is a non-negligible fraction of
        # the precession period.
        sim_rk4 = sp.LlgSimulator(yig)
        sim_rk4.set_magnetization(1.0, 0.0, 0.0)
        sim_rk4.set_external_field(0.0, 0.0, 0.1)
        traj_rk4 = sim_rk4.evolve(1e-9, 100, method="rk4")

        sim_euler = sp.LlgSimulator(yig)
        sim_euler.set_magnetization(1.0, 0.0, 0.0)
        sim_euler.set_external_field(0.0, 0.0, 0.1)
        traj_euler = sim_euler.evolve(1e-9, 100, method="euler")

        _, mx_rk4, _, _ = traj_rk4[-1]
        _, mx_euler, _, _ = traj_euler[-1]
        assert abs(mx_rk4 - mx_euler) > 0.05

    def test_norm_conserved_along_rk4_trajectory(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_magnetization(1.0, 0.0, 0.0)
        sim.set_external_field(0.0, 0.0, 0.1)
        traj = sim.evolve(1e-9, 500, method="rk4")
        for _, mx, my, mz in traj:
            norm = math.sqrt(mx * mx + my * my + mz * mz)
            assert norm == pytest.approx(1.0, abs=1e-9)

    def test_norm_conserved_along_euler_trajectory(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_magnetization(1.0, 0.0, 0.0)
        sim.set_external_field(0.0, 0.0, 0.1)
        traj = sim.evolve(1e-9, 500, method="euler")
        for _, mx, my, mz in traj:
            norm = math.sqrt(mx * mx + my * my + mz * mz)
            assert norm == pytest.approx(1.0, abs=1e-9)

    def test_evolve_with_zero_field_keeps_magnetization_constant(self, yig):
        sim = sp.LlgSimulator(yig)
        sim.set_magnetization(1.0, 0.0, 0.0)
        sim.set_external_field(0.0, 0.0, 0.0)
        traj = sim.evolve(1e-9, 50)
        for _, mx, my, mz in traj:
            assert (mx, my, mz) == pytest.approx((1.0, 0.0, 0.0), abs=1e-12)
