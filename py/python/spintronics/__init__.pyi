"""
Type stubs for the spintronics Python extension module.

These stubs cover all classes and constants exposed by the Rust/PyO3 bindings.
"""

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
