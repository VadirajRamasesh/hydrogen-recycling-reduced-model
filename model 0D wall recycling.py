"""
WEST / ITER reduced-order wall recycling model
Vadiraj, M.Sc. Hydrogen Technology, TH Rosenheim
December 2025

Reduced-order 0D reservoir model for wall inventory evolution.
Physics parameters (Wall Temp, Pumping) calibrated to WEST steady-state 
conditions [Loarer et al., Nucl. Fusion 60 126046 (2020)].

NOTE: This simulation applies a scaled particle flux to predict the 
theoretical breakdown of density control (R_eff > 1) in extended 
high-power discharges (t > 100s).
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass

k_B = 8.617333262145e-5

@dataclass(frozen=True)
class p:
    # Calibrated to WEST Physics (Loarer et al. 2020) 
    [cite_start]T0     = 573.0        # Max surface temp ~300°C (Actively cooled W) [cite: 171, 245]
    [cite_start]S_pump = 1e18         # Negligible active pumping in this config (~0.01 Pam^3/s) [cite: 257]
    
    # Tuned for Saturation Stress Test 
    R      = 0.992        # Global recycling coefficient (High recycling regime)
    S_in   = 1.2e22       # SCALED FLUX:Increased to force rapid saturation
    Nw_max = 1.15e23      # Effective wall capacity (Tuned to W coating volume)
    
    # Physics Constants
    tau_p  = 3.5          # Characteristic particle confinement time
    tau_0  = 1e-12        # Lattice vibration frequency factor
    E_a    = 1.05         # Effective de-trapping energy (eV)

def get_T_wall(t):
    if t < 50: return p.T0
    # Cubic approximation for surface heating excursion during long pulse
    # (Matches qualitative thermal rise leading to outgassing bursts)
    return p.T0 + 820 * ((t-50)/45)**3.1

def get_tau_w(t):
    # Arrhenius release time: tau = tau_0 * exp(Ea / kT)
    return p.tau_0 * np.exp(p.E_a / (k_B * get_T_wall(t)))

def rhs(t, y):
    Np, Nw = y
    if Np < 5e19: Np = 5e19 # Floor for numerical stability

    Gamma    = Np / p.tau_p
    prompt   = p.R * Gamma
    entering = (1-p.R) * Gamma
    
    # Wall uptake saturates as inventory Nw approaches Nw_max
    uptake   = max(0.0, 1.0 - Nw/p.Nw_max)
    
    # Thermal release driven by wall temperature evolution
    release  = Nw / get_tau_w(t)
    
    # Gas puff fueling profile (tanh ramp-up)
    fueling  = p.S_in * (1.0 + 1.8 * (1.0 + np.tanh((t-15)/4.0))/2.0)

    # 0D Balance Equations:
    # dNp/dt = Fueling + Recycled(Prompt + Thermal) -Loss(Gamma) - Pumping
    dNp_dt = fueling + prompt + release - Gamma - p.S_pump
    
    # dNw/dt = Entering_Flux * Uptake_Factor - Thermal_Release
    dNw_dt = entering * uptake - release
    
    return [dNp_dt, dNw_dt]

# Solve for 180s discharge to see crash
sol = solve_ivp(rhs, (0, 180), [1.8e21, 3.2e22], method='Radau',
                rtol=1e-8, atol=1e12, max_step=1.5)

# Post-processing calculations
Gamma_flux   = sol.y[0] / p.tau_p
release_flux = sol.y[1] / np.array([get_tau_w(ti) for ti in sol.t])

# Effective Recycling = (Prompt + Thermal) / Incident
R_eff        = (p.R * Gamma_flux + release_flux) / (Gamma_flux + 1e16)

# Plotting 
plt.figure(figsize=(10,11))

plt.subplot(411)
plt.plot(sol.t, sol.y[0]/1e21, 'tab:blue', lw=2.2)
plt.ylabel('$N_p$ [$10^{21}$]')
plt.title('Predicted Density Control Loss (WEST-Relevant Parameters)')
plt.grid(alpha=0.3)

plt.subplot(412)
plt.plot(sol.t, sol.y[1]/1e22, 'tab:orange', lw=2.2)
plt.axhline(p.Nw_max/1e22, color='r', ls='--', lw=1.8, label='Saturation Limit')
plt.ylabel('$N_w$ [$10^{22}$]')
plt.legend(loc='lower right')
plt.grid(alpha=0.3)

plt.subplot(413)
plt.plot(sol.t, (p.R * Gamma_flux + release_flux)/1e21, 'tab:green', lw=2.2)
plt.ylabel('Return flux [$10^{21}$ s$^{-1}$]')
plt.grid(alpha=0.3)

plt.subplot(414)
plt.plot(sol.t, R_eff, 'darkred', lw=2.8)
plt.axhline(0.995, color='k', ls='--')
plt.text(102, 0.996, 'LOSS OF DENSITY CONTROL', color='darkred', fontsize=12, weight='bold')
plt.ylabel('$R_{eff}$')
plt.xlabel('Time [s]')
plt.ylim(0.985, 1.002)
plt.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('WEST_Recycling_Model_Output.png', dpi=300, bbox_inches='tight')
plt.show()

