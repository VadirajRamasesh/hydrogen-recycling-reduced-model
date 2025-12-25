WEST / ITER reduced-order wall recycling model
Vadiraj, M.Sc. Hydrogen Technology, TH Rosenheim
December 2025

Reduced Model I tried to create  during Christmas break 2025

still needs proper multi-zone but good enough for now

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass

k_B = 8.617333262145e-5

@dataclass(frozen=True)
class p:
    tau_p  = 3.5          # WEST high-power phase
    tau_0  = 1e-12
    E_a    = 1.05         # Representative activation energy consistent with WEST tungsten studies (order ~1 eV)
    R      = 0.992          # effective baseline recycling coefficient (assumed/tuned for this scenario)
    Nw_max = 1.15e23
    S_in   = 1.2e22
    S_pump = 9e20
    T0     = 590.0

def get_T_wall(t):
    if t < 50: return p.T0
    # Simplified temperature ramp to reproduce qualitative thermal-release behavior (demonstration forcing).
    return p.T0 + 820 * ((t-50)/45)**3.1

def get_tau_w(t):
    return p.tau_0 * np.exp(p.E_a / (k_B * get_T_wall(t)))

def rhs(t, y):
    Np, Nw = y
    if Np < 5e19: Np = 5e19

    Gamma    = Np / p.tau_p
    prompt   = p.R * Gamma
    entering = (1-p.R) * Gamma
    uptake   = max(0.0, 1.0 - Nw/p.Nw_max)
    release  = Nw / get_tau_w(t)
    fueling  = p.S_in * (1.0 + 1.8 * (1.0 + np.tanh((t-15)/4.0))/2.0)

    return [fueling + prompt + release - Gamma - p.S_pump,
            entering * uptake - release]

sol = solve_ivp(rhs, (0, 180), [1.8e21, 3.2e22], method='Radau',
                rtol=1e-8, atol=1e12, max_step=1.5)

Gamma_flux   = sol.y[0] / p.tau_p
release_flux = sol.y[1] / np.array([get_tau_w(ti) for ti in sol.t])
R_eff        = (p.R * Gamma_flux + release_flux) / (Gamma_flux + 1e16)

plt.figure(figsize=(10,11))

plt.subplot(411)
plt.plot(sol.t, sol.y[0]/1e21, 'tab:blue', lw=2.2)
plt.ylabel('$N_p$ [$10^{21}$]')
plt.grid(alpha=0.3)

plt.subplot(412)
plt.plot(sol.t, sol.y[1]/1e22, 'tab:orange', lw=2.2)
plt.axhline(p.Nw_max/1e22, color='r', ls='--', lw=1.8)
plt.ylabel('$N_w$ [$10^{22}$]')
plt.grid(alpha=0.3)

plt.subplot(413)
plt.plot(sol.t, (p.R * Gamma_flux + release_flux)/1e21, 'tab:green', lw=2.2)
plt.ylabel('Return flux [$10^{21}$ s$^{-1}$]')
plt.grid(alpha=0.3)

plt.subplot(414)
plt.plot(sol.t, R_eff, 'darkred', lw=2.8)
plt.axhline(0.995, color='k', ls='--')
plt.text(102, 0.996, 'R_eff → 1 (loss of control leverage)', color='darkred', fontsize=12, weight='bold')
plt.ylabel('$R_{eff}$')
plt.xlabel('Time [s]')
plt.ylim(0.985, 1.002)
plt.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('Vadiraj_WEST_ITER_2025.png', dpi=300, bbox_inches='tight')
plt.show()

