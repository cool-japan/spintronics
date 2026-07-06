"""
spintronics - High-performance spintronics simulation library

Python bindings for the Rust spintronics library providing:
- Vector3: 3D vector operations for spintronics calculations
- LLG dynamics simulation (Landau-Lifshitz-Gilbert equation)
- LLB finite-temperature dynamics (Landau-Lifshitz-Bloch equation)
- Spin caloritronics (Onsager transport matrix, heat currents)
- Spin wave dispersion and magnon physics
- SIMD batch evolution for large spin ensembles
- Stochastic LLG with thermal noise
- Material presets (YIG, Permalloy, CoFeB, Iron, Cobalt, Nickel, ...)

Example usage::

    import spintronics

    # Create a ferromagnetic material
    yig = spintronics.Ferromagnet.yig()

    # Set up an LLG simulation
    sim = spintronics.LlgSimulator(yig)
    sim.set_magnetization(1.0, 0.0, 0.0)
    sim.set_external_field(0.0, 0.0, 0.1)  # 100 mT along z

    # Evolve for 1 ns
    trajectory = sim.evolve(1e-9, 1000)
    print(f"Final magnetization: {trajectory[-1]}")

    # Spin pumping simulation
    sp = spintronics.SpinPumpingSimulation()
    sp.set_fmr_conditions(9.65e9, 0.1)
    result = sp.run(1e-9, 1000)
    print(f"Peak ISHE voltage: {result['peak_voltage']:.4e} V")
"""
from .spintronics import *  # noqa: F401, F403
from .spintronics import (
    Vector3,
    Ferromagnet,
    SpinInterface,
    InverseSpinHall,
    LlgSimulator,
    SpinPumpingSimulation,
    LlbMaterial,
    LlbSolver,
    OnsagerMatrix,
    SpinCaloritronicsMaterial,
    batch_rk4_step,
    batch_rk4_multistep,
    HBAR,
    GAMMA,
    E_CHARGE,
    MU_B,
    KB,
)

__version__ = "0.3.2"
__all__ = [
    "Vector3",
    "Ferromagnet",
    "SpinInterface",
    "InverseSpinHall",
    "LlgSimulator",
    "SpinPumpingSimulation",
    "LlbMaterial",
    "LlbSolver",
    "OnsagerMatrix",
    "SpinCaloritronicsMaterial",
    "batch_rk4_step",
    "batch_rk4_multistep",
    "HBAR",
    "GAMMA",
    "E_CHARGE",
    "MU_B",
    "KB",
]
