# spintronics

Python bindings for the [spintronics](https://github.com/cool-japan/spintronics) Rust library — a pure-Rust, high-performance library for simulating spin dynamics, spin current generation, and conversion phenomena in magnetic and topological materials.

## Installation

```bash
pip install spintronics
```

## Quick Start

```python
import spintronics

# Create a ferromagnetic material (YIG preset)
yig = spintronics.Ferromagnet.yig()
print(f"YIG saturation magnetization: {yig.ms:.2e} A/m")
print(f"YIG Gilbert damping: {yig.alpha}")

# 3D vector operations
v1 = spintronics.Vector3(1.0, 0.0, 0.0)
v2 = spintronics.Vector3(0.0, 1.0, 0.0)
cross = v1.cross(v2)
print(f"Cross product: {cross}")

# LLG magnetization dynamics
sim = spintronics.LlgSimulator(yig)
sim.set_magnetization(1.0, 0.0, 0.0)
sim.set_external_field(0.0, 0.0, 0.1)  # 100 mT along z
trajectory = sim.evolve(1e-9, 1000)
print(f"Final state: t={trajectory[-1][0]:.2e} s, m=({trajectory[-1][1]:.4f}, {trajectory[-1][2]:.4f}, {trajectory[-1][3]:.4f})")

# Spin pumping + ISHE detection
sp = spintronics.SpinPumpingSimulation()
sp.set_fmr_conditions(9.65e9, 0.1)
result = sp.run(1e-9, 1000)
print(f"Peak ISHE voltage: {result['peak_voltage']:.4e} V")

# Physical constants
print(f"HBAR = {spintronics.HBAR:.4e} J·s")
print(f"GAMMA = {spintronics.GAMMA:.4e} rad/(s·T)")
print(f"KB = {spintronics.KB:.4e} J/K")
```

## Available Classes

| Class | Description |
|-------|-------------|
| `Vector3` | 3D vector with cross/dot product, normalization, arithmetic |
| `Ferromagnet` | Material parameters (damping, Ms, anisotropy); presets: YIG, Permalloy, CoFe, CoFeB, Fe, Co, Ni |
| `SpinInterface` | FM/NM interface (spin mixing conductance); presets: YIG/Pt, Py/Pt |
| `InverseSpinHall` | ISHE spin-to-charge conversion; presets: Pt, Ta, W |
| `LlgSimulator` | LLG equation integrator (RK4/Euler); trajectory output |
| `SpinPumpingSimulation` | Complete FMR → spin current → ISHE voltage workflow |

## Physical Constants

| Constant | Description |
|----------|-------------|
| `HBAR` | Reduced Planck constant (J·s) |
| `GAMMA` | Electron gyromagnetic ratio (rad/(s·T)) |
| `E_CHARGE` | Elementary charge (C) |
| `MU_B` | Bohr magneton (J/T) |
| `KB` | Boltzmann constant (J/K) |

## Building from Source

Requires Rust and [maturin](https://github.com/PyO3/maturin).

```bash
git clone https://github.com/cool-japan/spintronics
cd spintronics/py
pip install maturin
maturin develop --release
```

## License

MIT OR Apache-2.0

Copyright (c) 2025 COOLJAPAN OÜ (Team KitaSan)
