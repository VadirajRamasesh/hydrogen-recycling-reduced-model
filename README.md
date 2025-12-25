# Reduced-Order Hydrogen Recycling Model (0D)

Minimal reduced-order (0D) model for hydrogen recycling and wall inventory evolution in long-pulse operation.

## What this is
A two-state ODE model:
- `Np(t)` — lumped plasma particle inventory
- `Nw(t)` — lumped wall inventory

Processes represented:
- external fueling
- prompt recycling (baseline coefficient `R0`)
- wall uptake with finite capacity (`Nw_max`)
- thermally activated wall release via an Arrhenius-type residence time
- effective exhaust/pumping sink

**Purpose:** qualitative exploration and sensitivity studies (assumptions, time scales, regime transitions).  
**Not a validated predictive model.**

## Interpretation
The model outputs an effective recycling diagnostic:

\[
R_\mathrm{eff} = \frac{\text{prompt return} + \text{thermal release}}{\text{incident flux}}
\]

Any threshold line (e.g., `R_eff = 0.995`) is a **heuristic marker** to visualize approach to near-unity recycling, not a validated operational limit.

## Relationship to WEST / ITER literature
The plasma inventory balance follows the spirit of global particle balance reasoning used in tokamak long-pulse interpretation (e.g., WEST studies). Wall retention/release is represented by a simplified phenomenological closure to keep the model transparent and lightweight.

## How to run
```bash
python model_0d_recycling.py

