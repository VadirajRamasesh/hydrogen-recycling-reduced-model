"""
WEST / ITER reduced-order recycling model (0D)
Parameters consistent with Loarer et al. (Nucl. Fusion 2020) steady-state conditions.
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass, replace

# Boltzmann constant (eV/K)
k_B = 8.617e-5

@dataclass
class Params:
    # WEST Reference Conditions (Loarer 2020)
    T0: float     = 573.0      # Wall temp ~300 C (Actively cooled W)
    S_pump: float = 1e18       # Negligible pumping (~0.01 Pam^3/s)
    
    # Stress Test Parameters
    S_in: float   = 1.2e22     # SCALED FLUX (100x exp) for saturation test
    R: float      = 0.992      # High recycling regime
    Nw_max: float = 1.15e23    # Effective wall capacity
    
    # Physics Constants
    tau_p: float  = 3.5
    tau_0: float  = 1e-12
    E_a: float    = 1.05       # Activation energy (eV)

# Nominal configuration
P_NOMINAL = Params()

def get_T_wall(t, p: Params):
    if t < 50: return p.T0
    # Cubic fit approximation of thermal excursion
    return p.T0 + 820 * ((t-50)/45)**3.1

def get_tau(t, p: Params):
    # Arrhenius release time: tau = tau0 * exp(Ea/kT)
    T = get_T_wall(t, p)
    return p.tau_0 * np.exp(p.E_a / (k_B * T))

def rhs(t, y, p: Params):
    Np, Nw = y
    # Floor to prevent negative density numerical errors
    if Np < 1e15: Np = 1e15
    
    Gamma    = Np / p.tau_p
    prompt   = p.R * Gamma
    entering = (1 - p.R) * Gamma
    
    # Saturation factor (0=full, 1=empty)
    uptake   = max(0.0, 1.0 - Nw/p.Nw_max)
    
    # Thermal release
    tau = get_tau(t, p)
    release  = Nw / tau
    
    # Gas puff (tanh ramp-up)
    fueling  = p.S_in * (1.0 + 1.8 * (1.0 + np.tanh((t-15)/4.0))/2.0)

    # 0D Balance Equations
    dNp = fueling + prompt + release - Gamma - p.S_pump
    dNw = entering * uptake - release
    
    return [dNp, dNw]

def run_nominal():
    print("Running Nominal Simulation...")
    p = P_NOMINAL
    
    sol = solve_ivp(rhs, (0, 180), [1.8e21, 3.2e22], 
                    args=(p,), method='Radau', rtol=1e-8)
    
    # Post-processing diagnostics
    time = sol.t
    T_w = np.array([get_T_wall(ti, p) for ti in time])
    tau_w = p.tau_0 * np.exp(p.E_a / (k_B * T_w))
    
    Gamma = sol.y[0] / p.tau_p
    release = sol.y[1] / tau_w
    R_eff = (p.R * Gamma + release) / (Gamma + 1e-16)
    
    # Plotting
    plt.figure(figsize=(10, 10))
    
    plt.subplot(311)
    plt.plot(time, sol.y[0]/1e21, label='$N_p$ (Plasma)', lw=2)
    plt.ylabel('Inventory [$10^{21}$]')
    plt.title('WEST-Relevant Wall Saturation Model')
    plt.grid(alpha=0.3)
    plt.legend()
    
    plt.subplot(312)
    plt.plot(time, sol.y[1]/1e22, color='tab:orange', label='$N_w$ (Wall)', lw=2)
    plt.axhline(p.Nw_max/1e22, color='r', ls='--', label='Max Capacity')
    plt.ylabel('Wall Inv. [$10^{22}$]')
    plt.grid(alpha=0.3)
    plt.legend()
    
    plt.subplot(313)
    plt.plot(time, R_eff, color='darkred', lw=2)
    plt.axhline(1.0, color='k', ls='--')
    plt.text(105, 1.002, "Loss of Density Control", color='darkred', weight='bold')
    plt.ylabel('$R_{eff}$')
    plt.xlabel('Time [s]')
    plt.ylim(0.98, 1.01)
    plt.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('WEST_Nominal.png', dpi=150)
    print("Saved 'WEST_Nominal.png'")

def run_sweep():
    print("Running Activation Energy (Ea) Sensitivity Sweep...")
    Ea_values = [0.95, 1.0, 1.05, 1.10, 1.15]
    
    plt.figure(figsize=(10, 6))
    
    for val in Ea_values:
        # Create isolated parameter set for this run
        p_run = replace(P_NOMINAL, E_a=val)
        
        sol = solve_ivp(rhs, (0, 180), [1.8e21, 3.2e22], 
                        args=(p_run,), method='Radau', rtol=1e-8)
        
        # Reconstruct R_eff
        time = sol.t
        Np, Nw = sol.y
        Gamma = Np / p_run.tau_p
        
        T_w = np.array([get_T_wall(ti, p_run) for ti in time])
        tau_plot = p_run.tau_0 * np.exp(val / (k_B * T_w))
        
        R_eff = (p_run.R * Gamma + (Nw/tau_plot)) / (Gamma + 1e-16)
        
        label = f"Ea = {val} eV"
        if val == 1.05: label += " (Nominal)"
        plt.plot(time, R_eff, lw=2, label=label)

    plt.axhline(1.0, color='k', ls='--', label="Unity Limit")
    plt.title("Sensitivity: Activation Energy vs Control Loss Time")
    plt.ylabel("$R_{eff}$")
    plt.xlabel("Time [s]")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.ylim(0.98, 1.02)
    plt.tight_layout()
    plt.savefig('WEST_Sensitivity.png', dpi=150)
    print("Saved 'WEST_Sensitivity.png'")

if __name__ == "__main__":
    run_nominal()
    run_sweep()
