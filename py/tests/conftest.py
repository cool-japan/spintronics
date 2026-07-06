"""Shared pytest fixtures for the spintronics Python binding test suite.

These fixtures wrap the most commonly used material/interface/effect
presets so individual test modules don't have to repeat the same
construction calls.
"""
import pytest

import spintronics as sp


@pytest.fixture
def yig():
    """YIG (Yttrium Iron Garnet): ultra-low-damping ferrimagnetic insulator.

    The canonical material for spin pumping and magnon transport
    experiments (alpha ~= 1e-4).
    """
    return sp.Ferromagnet.yig()


@pytest.fixture
def permalloy():
    """Permalloy (Ni80Fe20): soft ferromagnetic alloy with near-zero
    magnetostriction, used as a higher-damping contrast to YIG."""
    return sp.Ferromagnet.permalloy()


@pytest.fixture
def yig_pt_interface():
    """Canonical YIG/Pt spin interface used in spin pumping experiments."""
    return sp.SpinInterface.yig_pt()


@pytest.fixture
def platinum():
    """Platinum inverse spin Hall effect (ISHE) converter (theta_SH ~ 0.08)."""
    return sp.InverseSpinHall.platinum()
