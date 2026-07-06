"""
Type stubs for the spintronics Python extension module.

These stubs cover all classes and constants exposed by the Rust/PyO3 bindings.
"""

import numpy as np
from numpy.typing import NDArray

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

HBAR: float
"""Reduced Planck constant (J·s)."""

GAMMA: float
"""Gyromagnetic ratio of the electron (rad/(s·T))."""

E_CHARGE: float
"""Elementary charge (C)."""

MU_B: float
"""Bohr magneton (J/T)."""

KB: float
"""Boltzmann constant (J/K)."""


# ---------------------------------------------------------------------------
# Vector3
# ---------------------------------------------------------------------------

class Vector3:
    """A 3D vector for spintronics calculations."""

    def __init__(self, x: float, y: float, z: float) -> None:
        """Create a new 3D vector.

        Args:
            x: X component.
            y: Y component.
            z: Z component.
        """

    @property
    def x(self) -> float:
        """X component."""

    @x.setter
    def x(self, value: float) -> None: ...

    @property
    def y(self) -> float:
        """Y component."""

    @y.setter
    def y(self, value: float) -> None: ...

    @property
    def z(self) -> float:
        """Z component."""

    @z.setter
    def z(self, value: float) -> None: ...

    def dot(self, other: Vector3) -> float:
        """Calculate the dot product with another vector."""

    def cross(self, other: Vector3) -> Vector3:
        """Calculate the cross product with another vector."""

    def magnitude(self) -> float:
        """Calculate the magnitude (Euclidean norm) of the vector."""

    def normalize(self) -> Vector3:
        """Return a unit vector in the same direction."""

    def to_tuple(self) -> tuple[float, float, float]:
        """Convert to a Python tuple (x, y, z)."""

    def to_list(self) -> list[float]:
        """Convert to a Python list [x, y, z]."""

    def __add__(self, other: Vector3) -> Vector3: ...
    def __sub__(self, other: Vector3) -> Vector3: ...
    def __mul__(self, scalar: float) -> Vector3: ...
    def __rmul__(self, scalar: float) -> Vector3: ...
    def __repr__(self) -> str: ...

    @staticmethod
    def unit_x() -> Vector3:
        """Create the unit vector along x (1, 0, 0)."""

    @staticmethod
    def unit_y() -> Vector3:
        """Create the unit vector along y (0, 1, 0)."""

    @staticmethod
    def unit_z() -> Vector3:
        """Create the unit vector along z (0, 0, 1)."""

    @staticmethod
    def zero() -> Vector3:
        """Create the zero vector (0, 0, 0)."""


# ---------------------------------------------------------------------------
# Ferromagnet
# ---------------------------------------------------------------------------

class Ferromagnet:
    """Ferromagnetic material properties.

    Holds the key parameters that define a ferromagnetic material:
    Gilbert damping, saturation magnetization, uniaxial anisotropy,
    easy axis direction, and exchange stiffness.
    """

    def __init__(
        self,
        alpha: float,
        ms: float,
        anisotropy_k: float = 0.0,
        easy_axis: tuple[float, float, float] = (0.0, 0.0, 1.0),
        exchange_a: float = 1e-11,
    ) -> None:
        """Create a ferromagnet with custom parameters.

        Args:
            alpha: Gilbert damping parameter (dimensionless).
            ms: Saturation magnetization (A/m).
            anisotropy_k: Uniaxial anisotropy constant (J/m^3).
            easy_axis: Easy axis direction (tuple of 3 floats, will be normalised).
            exchange_a: Exchange stiffness (J/m).
        """

    @property
    def alpha(self) -> float:
        """Gilbert damping parameter (dimensionless)."""

    @property
    def ms(self) -> float:
        """Saturation magnetization (A/m)."""

    @property
    def anisotropy_k(self) -> float:
        """Uniaxial anisotropy constant (J/m^3)."""

    @property
    def easy_axis(self) -> Vector3:
        """Easy axis direction (normalised Vector3)."""

    @property
    def exchange_a(self) -> float:
        """Exchange stiffness (J/m)."""

    @staticmethod
    def yig() -> Ferromagnet:
        """Create YIG (Yttrium Iron Garnet) material.

        YIG is a ferrimagnetic insulator with very low damping (~2e-4),
        ideal for spin pumping and magnon transport.
        """

    @staticmethod
    def permalloy() -> Ferromagnet:
        """Create Permalloy (Ni80Fe20) material.

        Soft magnetic alloy with near-zero magnetostriction and high permeability.
        """

    @staticmethod
    def cofe() -> Ferromagnet:
        """Create CoFe (Cobalt-Iron) alloy.

        High saturation magnetisation, commonly used in spin-transfer torque devices.
        """

    @staticmethod
    def cofeb() -> Ferromagnet:
        """Create CoFeB (Cobalt-Iron-Boron) alloy.

        Widely used in magnetic tunnel junctions due to perpendicular magnetic anisotropy.
        """

    @staticmethod
    def iron() -> Ferromagnet:
        """Create Iron (Fe) material."""

    @staticmethod
    def cobalt() -> Ferromagnet:
        """Create Cobalt (Co) material."""

    @staticmethod
    def nickel() -> Ferromagnet:
        """Create Nickel (Ni) material."""

    def __repr__(self) -> str: ...


# ---------------------------------------------------------------------------
# SpinInterface
# ---------------------------------------------------------------------------

class SpinInterface:
    """Spin interface between a ferromagnet and a normal metal.

    Characterised by the spin mixing conductance, interface normal direction,
    and interface area.
    """

    def __init__(
        self,
        g_r: float,
        g_i: float = 0.0,
        normal: tuple[float, float, float] = (0.0, 1.0, 0.0),
        area: float = 1e-12,
    ) -> None:
        """Create a spin interface.

        Args:
            g_r: Real part of spin mixing conductance (1/(Ohm·m^2)).
            g_i: Imaginary part of spin mixing conductance (1/(Ohm·m^2)).
            normal: Interface normal direction (tuple of 3 floats).
            area: Interface area (m^2).
        """

    @property
    def g_r(self) -> float:
        """Real part of spin mixing conductance (1/(Ohm·m^2))."""

    @property
    def g_i(self) -> float:
        """Imaginary part of spin mixing conductance (1/(Ohm·m^2))."""

    @property
    def normal(self) -> Vector3:
        """Interface normal direction (Vector3)."""

    @property
    def area(self) -> float:
        """Interface area (m^2)."""

    @staticmethod
    def yig_pt() -> SpinInterface:
        """Create a YIG/Pt interface (canonical spin pumping system)."""

    @staticmethod
    def py_pt() -> SpinInterface:
        """Create a Permalloy/Pt interface."""

    def __repr__(self) -> str: ...


# ---------------------------------------------------------------------------
# InverseSpinHall
# ---------------------------------------------------------------------------

class InverseSpinHall:
    """Inverse Spin Hall Effect (ISHE) converter.

    Converts spin current to charge current via spin-orbit coupling.
    The generated electric field satisfies:
        E = rho * theta_SH * (j_s x sigma)
    """

    def __init__(self, theta_sh: float, rho: float) -> None:
        """Create an ISHE converter.

        Args:
            theta_sh: Spin Hall angle (dimensionless).
            rho: Electrical resistivity (Ohm·m).
        """

    @property
    def theta_sh(self) -> float:
        """Spin Hall angle (dimensionless)."""

    @property
    def rho(self) -> float:
        """Electrical resistivity (Ohm·m)."""

    def convert(self, js_flow: Vector3, js_polarization: Vector3) -> Vector3:
        """Convert spin current to an electric field vector.

        Args:
            js_flow: Spin current flow direction (Vector3).
            js_polarization: Spin polarisation with magnitude (Vector3, A/m^2).

        Returns:
            Electric field vector (V/m).
        """

    def voltage(
        self,
        js_flow: Vector3,
        js_polarization: Vector3,
        strip_width: float,
    ) -> float:
        """Calculate the ISHE voltage across a strip.

        Args:
            js_flow: Spin current flow direction (Vector3).
            js_polarization: Spin polarisation with magnitude (Vector3, A/m^2).
            strip_width: Width of the conducting strip (m).

        Returns:
            Voltage (V).
        """

    def efficiency(self) -> float:
        """Return the conversion efficiency (V/W ratio)."""

    @staticmethod
    def platinum() -> InverseSpinHall:
        """Create a Platinum (Pt) ISHE converter (theta_SH ~ +0.08)."""

    @staticmethod
    def tantalum() -> InverseSpinHall:
        """Create a Tantalum (Ta) ISHE converter (theta_SH ~ +0.12)."""

    @staticmethod
    def tungsten() -> InverseSpinHall:
        """Create a Tungsten (W) ISHE converter (theta_SH ~ -0.30)."""

    def __repr__(self) -> str: ...


# ---------------------------------------------------------------------------
# LlgSimulator
# ---------------------------------------------------------------------------

class LlgSimulator:
    """LLG (Landau-Lifshitz-Gilbert) equation simulator.

    Simulates single-macrospin magnetisation dynamics under an effective field:
        dm/dt = -gamma * (m x H_eff) + alpha * (m x dm/dt)
    """

    def __init__(self, material: Ferromagnet) -> None:
        """Create a new LLG simulator.

        Args:
            material: Ferromagnetic material parameters.
        """

    @property
    def time(self) -> float:
        """Current simulation time (seconds)."""

    def set_magnetization(self, mx: float, my: float, mz: float) -> None:
        """Set the magnetisation direction (automatically normalised).

        Args:
            mx, my, mz: Magnetisation components.
        """

    def get_magnetization(self) -> Vector3:
        """Return the current magnetisation as a Vector3."""

    def set_external_field(self, hx: float, hy: float, hz: float) -> None:
        """Set the external magnetic field.

        Args:
            hx, hy, hz: Field components (Tesla).
        """

    def get_external_field(self) -> Vector3:
        """Return the current external field as a Vector3."""

    def reset_time(self) -> None:
        """Reset simulation time to zero."""

    def dm_dt(self) -> Vector3:
        """Calculate dm/dt at the current state."""

    def step_rk4(self, dt: float) -> None:
        """Perform one RK4 integration step.

        Args:
            dt: Time step (seconds).
        """

    def step_euler(self, dt: float) -> None:
        """Perform one Euler integration step.

        Args:
            dt: Time step (seconds).
        """

    def evolve(
        self,
        duration: float,
        n_steps: int,
        method: str = "rk4",
    ) -> list[tuple[float, float, float, float]]:
        """Evolve magnetisation for a specified duration.

        Args:
            duration: Total simulation time (seconds).
            n_steps: Number of integration steps.
            method: Integration method — ``"rk4"`` (default) or ``"euler"``.

        Returns:
            List of ``(time, mx, my, mz)`` tuples, one per recorded step.
        """

    def precession_frequency(self) -> float:
        """Return the Larmor precession frequency (rad/s)."""

    def precession_period(self) -> float:
        """Return the precession period (seconds)."""

    def __repr__(self) -> str: ...


# ---------------------------------------------------------------------------
# SpinPumpingSimulation
# ---------------------------------------------------------------------------

class SpinPumpingSimulation:
    """Complete spin pumping simulation (FMR -> spin current -> ISHE voltage).

    Implements the canonical YIG/Pt experiment workflow:
    1. Magnetisation precession under FMR driving.
    2. Spin pumping current generation at the FM/NM interface.
    3. ISHE voltage detection in the normal metal layer.
    """

    def __init__(self) -> None:
        """Create a spin pumping simulation with default YIG/Pt parameters."""

    def set_sample_length(self, length: float) -> None:
        """Set the sample length for voltage measurement (m)."""

    def set_field(self, hx: float, hy: float, hz: float) -> None:
        """Set the external magnetic field (Tesla)."""

    def set_magnetization(self, mx: float, my: float, mz: float) -> None:
        """Set the initial magnetisation direction (automatically normalised)."""

    def set_fmr_conditions(self, frequency: float, field: float) -> None:
        """Configure FMR driving conditions.

        Args:
            frequency: Microwave frequency (Hz).
            field: DC magnetic field magnitude (T), applied along z.
        """

    def run(self, duration: float, n_steps: int) -> dict:
        """Run the simulation.

        Args:
            duration: Total simulation time (seconds).
            n_steps: Number of integration steps.

        Returns:
            Dictionary containing:
            - ``times``: list of time values (s)
            - ``mx``, ``my``, ``mz``: magnetisation components
            - ``spin_current``: spin current magnitude at each step (A/m^2)
            - ``voltage``: ISHE voltage at each step (V)
            - ``peak_voltage``: maximum ISHE voltage (V)
            - ``avg_voltage``: time-averaged ISHE voltage (V)
            - ``peak_spin_current``: maximum spin current magnitude (A/m^2)
        """

    def get_magnetization(self) -> Vector3:
        """Return the current magnetisation as a Vector3."""

    def __repr__(self) -> str: ...


# ---------------------------------------------------------------------------
# LlbMaterial
# ---------------------------------------------------------------------------

class LlbMaterial:
    """Material parameters for the LLB (Landau-Lifshitz-Bloch) equation.

    Encapsulates Curie temperature, Gilbert damping, quantum spin number,
    and zero-temperature saturation magnetization needed for finite-temperature
    spin dynamics.
    """

    def __init__(
        self,
        curie_temp: float,
        alpha: float,
        spin_s: float,
        ms_0: float,
    ) -> None:
        """Create a new LLB material with explicit parameters.

        Args:
            curie_temp: Curie temperature T_C (K).
            alpha: Base Gilbert damping constant (dimensionless).
            spin_s: Quantum spin number S (e.g. 0.5 for spin-1/2).
            ms_0: Zero-temperature saturation magnetization (A/m).
        """

    @property
    def curie_temp(self) -> float:
        """Curie temperature T_C (K)."""

    @property
    def alpha(self) -> float:
        """Base Gilbert damping constant alpha (dimensionless)."""

    @property
    def spin_s(self) -> float:
        """Quantum spin number S."""

    @property
    def ms_0(self) -> float:
        """Zero-temperature saturation magnetization M_s (A/m)."""

    def equilibrium_magnetization(self, temperature: float) -> float:
        """Compute equilibrium magnetization m_e(T) via the self-consistent Brillouin equation.

        Returns 0.0 for T >= T_C or T <= 0.

        Args:
            temperature: Temperature (K).

        Returns:
            Dimensionless equilibrium magnetization in [0, 1].
        """

    def alpha_parallel(self, temperature: float) -> float:
        """Temperature-dependent longitudinal damping alpha_parallel(T).

        - T < T_C: alpha_parallel = alpha * (2/5 + 3*T/(5*T_C))
        - T >= T_C: alpha_parallel = 2*alpha * T/(5*T_C)

        Args:
            temperature: Temperature (K).
        """

    def alpha_perp(self, temperature: float) -> float:
        """Temperature-dependent transverse damping alpha_perp(T).

        - T < T_C: alpha_perp = alpha * T/T_C (minimum 1e-10)
        - T >= T_C: alpha_perp = 2*alpha * T/(5*T_C)

        Args:
            temperature: Temperature (K).
        """

    @staticmethod
    def iron() -> LlbMaterial:
        """Iron (Fe) LLB material preset (T_C = 1043 K, alpha = 0.01, S = 1.0)."""

    @staticmethod
    def nickel() -> LlbMaterial:
        """Nickel (Ni) LLB material preset (T_C = 631 K, alpha = 0.064, S = 0.3)."""

    @staticmethod
    def cofeb() -> LlbMaterial:
        """CoFeB LLB material preset (T_C = 1000 K, alpha = 0.005, S = 0.5)."""

    def __repr__(self) -> str: ...


# ---------------------------------------------------------------------------
# LlbSolver
# ---------------------------------------------------------------------------

class LlbSolver:
    """LLB equation solver with 4th-order Runge-Kutta integration.

    Unlike the standard LLG solver, the LLB solver allows |m| to change
    over time via longitudinal relaxation. This is essential near and above
    the Curie temperature T_C.
    """

    def __init__(
        self,
        material: LlbMaterial,
        dt: float,
        temperature: float,
        h_ext: tuple[float, float, float],
    ) -> None:
        """Create a new LLB solver.

        Args:
            material: LlbMaterial with Curie temperature, damping, spin number, M_s.
            dt: Integration time step (s).
            temperature: Simulation temperature (K).
            h_ext: External applied field (hx, hy, hz) (T).
        """

    @property
    def dt(self) -> float:
        """Integration time step dt (s)."""

    @property
    def temperature(self) -> float:
        """Simulation temperature (K)."""

    @property
    def h_ext(self) -> tuple[float, float, float]:
        """External field h_ext as (hx, hy, hz) tuple (T)."""

    @property
    def gamma(self) -> float:
        """Gyromagnetic ratio (rad/(s·T))."""

    def step(self, m: tuple[float, float, float]) -> list[float]:
        """Advance the magnetization by one RK4 time step.

        Note: Unlike LLG, |m| is NOT renormalized after each step. The LLB
        equation explicitly evolves the magnetization magnitude.

        Args:
            m: Current magnetization vector (mx, my, mz).

        Returns:
            Updated magnetization vector [mx, my, mz] after one dt.
        """

    def run(
        self,
        m0: tuple[float, float, float],
        num_steps: int,
        record_every: int,
    ) -> dict:
        """Run the LLB simulation for num_steps time steps.

        Args:
            m0: Initial magnetization vector (mx, my, mz).
            num_steps: Total number of integration steps.
            record_every: Snapshot interval (1 = every step, N = every Nth step).

        Returns:
            Dictionary containing:
            - ``mx``, ``my``, ``mz``: list of float — trajectory components
            - ``m_magnitude``: list of float — |m(t)| at each snapshot
            - ``time``: list of float — time stamps (s)
            - ``equilibrium_m``: float — m_e(T) for the simulation temperature
        """

    def __repr__(self) -> str: ...


# ---------------------------------------------------------------------------
# OnsagerMatrix
# ---------------------------------------------------------------------------

class OnsagerMatrix:
    """Onsager transport matrix for spin caloritronics.

    Encodes the linear-response coupling between charge, spin, and heat
    currents in a magnetic heterostructure using Onsager reciprocal relations.
    """

    def __init__(
        self,
        temperature: float,
        conductivity: float,
        seebeck: float,
        spin_seebeck: float,
        hall_angle: float,
        thermal_conductivity: float,
    ) -> None:
        """Create an Onsager matrix with explicit parameters.

        Args:
            temperature: Temperature (K).
            conductivity: Electrical conductivity sigma (S/m).
            seebeck: Seebeck coefficient S_e (V/K).
            spin_seebeck: Spin Seebeck coefficient S_s (A/(m·K)).
            hall_angle: Anomalous Hall angle theta_H (dimensionless).
            thermal_conductivity: Thermal conductivity kappa (W/(m·K)).
        """

    @property
    def temperature(self) -> float:
        """Temperature (K)."""

    @property
    def conductivity(self) -> float:
        """Electrical conductivity sigma (S/m)."""

    @property
    def seebeck(self) -> float:
        """Seebeck coefficient S_e (V/K)."""

    @property
    def spin_seebeck(self) -> float:
        """Spin Seebeck coefficient S_s (A/(m·K))."""

    @property
    def hall_angle(self) -> float:
        """Anomalous Hall angle theta_H (dimensionless)."""

    @property
    def thermal_conductivity(self) -> float:
        """Thermal conductivity kappa (W/(m·K))."""

    def reciprocity_error(self) -> float:
        """Check Onsager reciprocity: returns the relative deviation from L_ij = T * L_ji.

        For analytically constructed matrices this is zero by construction.
        Values < 1e-10 indicate reciprocity is satisfied.
        """

    def spin_current_from_grad_t(
        self, grad_t: tuple[float, float, float]
    ) -> list[float]:
        """Compute the spin current density from a temperature gradient (spin Seebeck effect).

        j_s = S_s * sigma * grad(T)  (A/m^2)

        Args:
            grad_t: Temperature gradient vector (dT/dx, dT/dy, dT/dz) (K/m).

        Returns:
            Spin current density [jx, jy, jz] (A/m^2).
        """

    def heat_current_from_spin_current(
        self, j_spin: tuple[float, float, float]
    ) -> list[float]:
        """Compute the heat current from a spin current (spin Peltier effect).

        j_Q = Pi_s * j_s = T * S_s * j_s  (W/m^2)

        Args:
            j_spin: Spin current density (jx, jy, jz) (A/m^2).

        Returns:
            Heat current density [jQx, jQy, jQz] (W/m^2).
        """

    def nernst_voltage(self, grad_t: float) -> float:
        """Compute the anomalous Nernst voltage per unit length.

        nu = -S_e * theta_H * |grad(T)|  (V/m)

        Args:
            grad_t: Longitudinal temperature gradient magnitude (K/m).

        Returns:
            Transverse Nernst voltage (V/m).
        """

    def all_currents(
        self,
        grad_t: tuple[float, float, float],
        e_field: tuple[float, float, float],
    ) -> dict:
        """Compute charge, spin, and heat currents simultaneously.

        Full Onsager response:
            j_c = sigma * (E + S_e * grad(T))      (A/m^2)
            j_s = S_s * sigma * grad(T)             (A/m^2)
            j_Q = T*S_e*sigma*E - kappa*grad(T)     (W/m^2)

        Args:
            grad_t: Temperature gradient (dT/dx, dT/dy, dT/dz) (K/m).
            e_field: Electric field (Ex, Ey, Ez) (V/m).

        Returns:
            Dictionary containing:
            - ``charge_current``: list of 3 floats — j_c (A/m^2)
            - ``spin_current``: list of 3 floats — j_s (A/m^2)
            - ``heat_current``: list of 3 floats — j_Q (W/m^2)
        """

    @staticmethod
    def yig_pt(temperature: float) -> OnsagerMatrix:
        """YIG/Pt bilayer preset at `temperature` (K)."""

    @staticmethod
    def fe_pt(temperature: float) -> OnsagerMatrix:
        """Fe/Pt bilayer preset at `temperature` (K)."""

    @staticmethod
    def cofeb_pt(temperature: float) -> OnsagerMatrix:
        """CoFeB/Pt bilayer preset at `temperature` (K)."""

    def __repr__(self) -> str: ...


# ---------------------------------------------------------------------------
# SpinCaloritronicsMaterial
# ---------------------------------------------------------------------------

class SpinCaloritronicsMaterial:
    """Unified spin caloritronic material: combines an Onsager matrix with a
    heat current calculator for a single high-level computation.

    This is the primary user-facing type for spin caloritronics calculations.
    Use ``compute_all()`` to obtain all cross-effects simultaneously.

    This class has no public constructor; obtain instances via the
    ``yig_pt``, ``fe_pt``, ``cofeb_pt``, or ``from_onsager`` factory methods.
    """

    @property
    def name(self) -> str:
        """Human-readable material system label (e.g. "YIG/Pt", "Fe/Pt")."""

    @property
    def onsager(self) -> OnsagerMatrix:
        """The underlying OnsagerMatrix as a Python object."""

    def compute_all(
        self,
        grad_t: tuple[float, float, float],
        j_spin: tuple[float, float, float],
    ) -> dict:
        """Compute all spin-caloritronic cross-effects.

        Evaluates:
        - Spin Seebeck current: j_s = S_s * sigma * grad(T)
        - Spin Peltier heat: |Pi_s * j_s|
        - Anomalous Nernst voltage: nu = -S_e * theta_H * |grad(T)|
        - Spin Nernst current (proxy): theta_H * j_s^SSE
        - Total heat current: -kappa*grad(T) + Pi_s*j_s
        - Onsager reciprocity check

        Args:
            grad_t: Temperature gradient (dT/dx, dT/dy, dT/dz) (K/m).
            j_spin: Injected spin current density (jx, jy, jz) (A/m^2).

        Returns:
            Dictionary containing:
            - ``spin_seebeck_current``: list of float — j_s (A/m^2)
            - ``peltier_heat``: float — |j_Q^sPeltier| (W/m^2)
            - ``nernst_voltage``: float — anomalous Nernst nu (V/m)
            - ``spin_nernst_current``: list of float — j_s^SN proxy (A/m^2)
            - ``reciprocity_satisfied``: bool — Onsager error < 1e-10
            - ``total_heat_current``: list of float — total j_Q (W/m^2)
        """

    @staticmethod
    def yig_pt(temperature: float) -> SpinCaloritronicsMaterial:
        """YIG/Pt bilayer preset at `temperature` (K)."""

    @staticmethod
    def fe_pt(temperature: float) -> SpinCaloritronicsMaterial:
        """Fe/Pt bilayer preset at `temperature` (K)."""

    @staticmethod
    def cofeb_pt(temperature: float) -> SpinCaloritronicsMaterial:
        """CoFeB/Pt bilayer preset at `temperature` (K)."""

    @staticmethod
    def from_onsager(onsager: OnsagerMatrix) -> SpinCaloritronicsMaterial:
        """Create a material from an Onsager matrix.

        The heat-current calculator is built automatically via Kelvin relations:
            Pi   = T * S_e
            Pi_s = T * S_s

        Args:
            onsager: OnsagerMatrix for the desired material system.
        """

    def __repr__(self) -> str: ...


# ---------------------------------------------------------------------------
# SIMD batch LLG evolution (numpy)
# ---------------------------------------------------------------------------

def batch_rk4_step(
    m: NDArray[np.float64],
    h_eff: NDArray[np.float64],
    alpha: float,
    gamma: float,
    dt: float,
) -> NDArray[np.float64]:
    """Evolve N spins for one RK4 time step using SIMD batch processing.

    Each spin evolves independently under its own effective field (no coupling).
    The result is normalized so that |m_i| = 1 for every spin.

    Args:
        m: (N, 3) numpy array of magnetization vectors (row-major, each row is m_i).
        h_eff: (N, 3) numpy array of effective field vectors (T).
        alpha: Gilbert damping constant (dimensionless).
        gamma: Gyromagnetic ratio (rad/(s·T)), typically about 1.761e11.
        dt: Integration time step (s).

    Returns:
        (N, 3) numpy array of updated magnetization vectors (unit length).

    Raises:
        ValueError: If array shapes are incompatible or lengths differ.
    """

def batch_rk4_multistep(
    m: NDArray[np.float64],
    h_eff: NDArray[np.float64],
    alpha: float,
    gamma: float,
    dt: float,
    num_steps: int,
) -> NDArray[np.float64]:
    """Evolve N spins for ``num_steps`` RK4 time steps using SIMD batch processing.

    ``h_eff`` is held constant throughout (static field approximation). This is
    more efficient than calling ``batch_rk4_step`` in a loop because no
    Python-to-Rust conversion overhead is incurred at each step.

    The result is normalized so that |m_i| = 1 after all steps.

    Args:
        m: (N, 3) numpy array of magnetization vectors (initial state).
        h_eff: (N, 3) numpy array of effective field vectors (T) (constant).
        alpha: Gilbert damping constant (dimensionless).
        gamma: Gyromagnetic ratio (rad/(s·T)).
        dt: Integration time step (s).
        num_steps: Number of RK4 iterations to perform.

    Returns:
        (N, 3) numpy array of final magnetization vectors (unit length).

    Raises:
        ValueError: If array shapes are incompatible or lengths differ.
    """
