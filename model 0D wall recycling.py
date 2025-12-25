"""
Reduced-order (0D) hydrogen recycling / wall inventory model.

Scope
- Minimal two-reservoir ODE system:
    Np(t): lumped plasma particle inventory   [particles]
    Nw(t): lumped wall inventory              [particles]
- Includes: fueling, prompt recycling, wall uptake with finite capacity,
  and thermally activated wall release via an Arrhenius-type residence time.
- Intended for qualitative behavior and sensitivity studies (not validation).

Notes
- Global particle-balance framing is consistent with approaches used to
  interpret long-pulse tokamak operation (e.g., WEST literature).
- Wall retention/release is represented by a simplified phenomenological closure.
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass

# Boltzmann constant in eV/K (use with Ea in eV and T in K)
k_B_eV = 8.617333262145e-5


@dataclass(frozen=True)
class Params:
    # All sources/sinks S_* are in [particles/s]. Inventories N* are [particles].

    # Scenario inputs (illustrative)
    T0_K: float = 573.0         # baseline wall temperature [K]
    S_in: float = 1.2e22        # fueling source term [particles/s]
    S_pump: float = 1.0e18      # effective exhaust sink [particles/s]

    # Recycling / wall inventory (illustrative)
    R0: float = 0.992           # baseline prompt recycling coefficient [-]
    Nw_max: float = 1.15e23     # effective wall capacity [particles]

    # Time scales (illustrative)
    tau_p: float = 3.5          # effective particle confinement time [s]

    # Thermally activated release (phenomenological closure)
    tau0: float = 1.0e-12       # pre-exponential time scale [s]
    Ea_eV: float = 1.05         # effective activation energy [eV]


P = Params()


def T_wall_K(t: float) -> float:
    """
    Prescribed wall temperature trajectory [K] (scenario forcing).
    Smooth ramp is used to trigger stronger thermally activated release later in time.
    """
    if t < 50.0:
        return P.T0_K
    return P.T0_K + 820.0 * ((t - 50.0) / 45.0) ** 3.1


def tau_release_s(t: float) -> float:
    """Arrhenius-type residence time: tau = tau0 * exp(Ea / (kB*T))."""
    return P.tau0 * np.exp(P.Ea_eV / (k_B_eV * T_wall_K(t)))


def rhs(t: float, y: np.ndarray) -> list[float]:
    Np, Nw = float(y[0]), float(y[1])

    # keep inventories non-negative
    Np = max(Np, 0.0)
    Nw = max(Nw, 0.0)

    # incident flux proxy
    Gamma = Np / P.tau_p  # [particles/s]

    # prompt recycling and wall-entering flux
    prompt = P.R0 * Gamma
    entering = (1.0 - P.R0) * Gamma

    # wall uptake saturates as Nw -> Nw_max
    uptake = np.clip(1.0 - Nw / P.Nw_max, 0.0, 1.0)

    # thermally activated wall release
    release = Nw / tau_release_s(t)

    # fueling profile (smooth ramp)
    fueling = P.S_in * (1.0 + 1.8 * (1.0 + np.tanh((t - 15.0) / 4.0)) / 2.0)

    dNp_dt = fueling + prompt + release - Gamma - P.S_pump
    dNw_dt = entering * uptake - release

    return [dNp_dt, dNw_dt]


def main() -> None:
    t0, tf = 0.0, 180.0
    y0 = [1.8e21, 3.2e22]  # initial inventories [particles]

    sol = solve_ivp(
        rhs, (t0, tf), y0,
        method="Radau",
        rtol=1e-8,
        atol=1e12,
        max_step=1.5
    )

    t = sol.t
    Np = sol.y[0]
    Nw = sol.y[1]

    Gamma = Np / P.tau_p
    tauw = np.array([tau_release_s(ti) for ti in t])
    release_flux = Nw / tauw

    # Effective recycling diagnostic: returned / incident (prompt + thermal) / incident
    eps = 1e-30
    R_eff = (P.R0 * Gamma + release_flux) / (Gamma + eps)

    # heuristic marker only (for visualization)
    R_mark = 0.995

    plt.figure(figsize=(10, 11))

    plt.subplot(411)
    plt.plot(t, Np / 1e21, lw=2.0)
    plt.ylabel(r"$N_p$ [$10^{21}$]")
    plt.grid(alpha=0.3)

    plt.subplot(412)
    plt.plot(t, Nw / 1e22, lw=2.0)
    plt.axhline(P.Nw_max / 1e22, ls="--", lw=1.5)
    plt.ylabel(r"$N_w$ [$10^{22}$]")
    plt.grid(alpha=0.3)

    plt.subplot(413)
    plt.plot(t, (P.R0 * Gamma + release_flux) / 1e21, lw=2.0)
    plt.ylabel(r"Return flux [$10^{21}$ s$^{-1}$]")
    plt.grid(alpha=0.3)

    plt.subplot(414)
    plt.plot(t, R_eff, lw=2.2)
    plt.axhline(R_mark, ls="--", lw=1.2)
    plt.text(102, R_mark + 0.0007, "near-unity recycling marker (heuristic)", fontsize=11)
    plt.ylabel(r"$R_{\mathrm{eff}}$")
    plt.xlabel("Time [s]")
    plt.ylim(0.985, 1.01)
    plt.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig("WEST_Recycling_Model_Output.png", dpi=300, bbox_inches="tight")
    plt.show()

    print("Run complete.")
    print(f"Final inventories: Np={Np[-1]:.3e}, Nw={Nw[-1]:.3e}")
    print(f"Final R_eff: {R_eff[-1]:.6f}")


if __name__ == "__main__":
    main()

