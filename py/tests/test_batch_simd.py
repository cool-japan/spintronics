"""Bonus tests for the SIMD batch LLG evolution functions: batch_rk4_step,
batch_rk4_multistep (src/python/batch.rs), which operate on (N, 3) numpy
arrays in a structure-of-arrays layout internally.

See the module docstring in test_llb.py for a note on why these
functions are tested here despite not yet appearing in `__all__` / the
`.pyi` stub.
"""
import numpy as np
import pytest

import spintronics as sp


def _make_batch(n, mx=1.0, hz=0.1):
    m = np.zeros((n, 3))
    m[:, 0] = mx
    h = np.zeros((n, 3))
    h[:, 2] = hz
    return m, h


class TestBatchRk4Step:
    def test_returns_expected_shape_and_dtype(self):
        m, h = _make_batch(8)
        result = sp.batch_rk4_step(m, h, alpha=0.01, gamma=sp.GAMMA, dt=1e-13)
        assert result.shape == (8, 3)
        assert result.dtype == np.float64

    def test_result_is_unit_normalized(self):
        m, h = _make_batch(16)
        result = sp.batch_rk4_step(m, h, alpha=0.01, gamma=sp.GAMMA, dt=1e-13)
        norms = np.linalg.norm(result, axis=1)
        assert norms == pytest.approx(np.ones(16))

    def test_matches_single_spin_llg_simulator_step(self):
        m, h = _make_batch(1)
        result = sp.batch_rk4_step(m, h, alpha=0.01, gamma=sp.GAMMA, dt=1e-13)

        fm = sp.Ferromagnet(alpha=0.01, ms=1e5)
        sim = sp.LlgSimulator(fm)
        sim.set_magnetization(1.0, 0.0, 0.0)
        sim.set_external_field(0.0, 0.0, 0.1)
        sim.step_rk4(1e-13)

        assert tuple(result[0]) == pytest.approx(
            sim.get_magnetization().to_tuple(), abs=1e-9
        )

    def test_mismatched_row_counts_raise_value_error(self):
        m, _ = _make_batch(8)
        _, h_bad = _make_batch(9)
        with pytest.raises(ValueError):
            sp.batch_rk4_step(m, h_bad, alpha=0.01, gamma=sp.GAMMA, dt=1e-13)

    def test_wrong_number_of_columns_raises_value_error(self):
        m, _ = _make_batch(8)
        bad_h = np.zeros((8, 4))
        with pytest.raises(ValueError):
            sp.batch_rk4_step(m, bad_h, alpha=0.01, gamma=sp.GAMMA, dt=1e-13)


class TestBatchRk4Multistep:
    def test_returns_expected_shape(self):
        m, h = _make_batch(8)
        result = sp.batch_rk4_multistep(
            m, h, alpha=0.01, gamma=sp.GAMMA, dt=1e-13, num_steps=50
        )
        assert result.shape == (8, 3)

    def test_result_is_unit_normalized(self):
        m, h = _make_batch(8)
        result = sp.batch_rk4_multistep(
            m, h, alpha=0.01, gamma=sp.GAMMA, dt=1e-13, num_steps=50
        )
        norms = np.linalg.norm(result, axis=1)
        assert norms == pytest.approx(np.ones(8))

    def test_matches_repeated_single_step_calls(self):
        m, h = _make_batch(4)
        num_steps = 20
        dt = 1e-13

        multistep_result = sp.batch_rk4_multistep(
            m, h, alpha=0.01, gamma=sp.GAMMA, dt=dt, num_steps=num_steps
        )

        looped = m.copy()
        for _ in range(num_steps):
            looped = sp.batch_rk4_step(looped, h, alpha=0.01, gamma=sp.GAMMA, dt=dt)

        assert multistep_result == pytest.approx(looped, abs=1e-9)

    def test_matches_llg_simulator_over_many_steps(self):
        m, h = _make_batch(1)
        num_steps = 100
        dt = 1e-13
        result = sp.batch_rk4_multistep(
            m, h, alpha=0.01, gamma=sp.GAMMA, dt=dt, num_steps=num_steps
        )

        fm = sp.Ferromagnet(alpha=0.01, ms=1e5)
        sim = sp.LlgSimulator(fm)
        sim.set_magnetization(1.0, 0.0, 0.0)
        sim.set_external_field(0.0, 0.0, 0.1)
        for _ in range(num_steps):
            sim.step_rk4(dt)

        assert tuple(result[0]) == pytest.approx(
            sim.get_magnetization().to_tuple(), abs=1e-9
        )

    def test_mismatched_row_counts_raise_value_error(self):
        m, _ = _make_batch(4)
        _, h_bad = _make_batch(5)
        with pytest.raises(ValueError):
            sp.batch_rk4_multistep(
                m, h_bad, alpha=0.01, gamma=sp.GAMMA, dt=1e-13, num_steps=10
            )
