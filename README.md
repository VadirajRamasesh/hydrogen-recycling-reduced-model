 Reduced-Order Wall Recycling Model for Hydrogen in Fusion Devices

 Overview
This repository contains a physics-based reduced-order (0D) model for hydrogen recycling at plasma-facing materials, developed to study global flux balance and wall inventory evolution in magnetic fusion devices.

The model represents the coupled dynamics of plasma particle inventory and wall hydrogen inventory using lumped balance equations, without spatial resolution inside the material.

 Physical Model
The model tracks two time-dependent state variables:
- Plasma hydrogen inventory
- Wall hydrogen inventory

The governing balances include:
- Prompt recycling at the plasma–wall interface
- Wall uptake limited by finite storage capacity
- Delayed hydrogen release from the wall
- External fueling and pumping terms

Wall release is modeled using a temperature-dependent effective residence time with Arrhenius behavior.

 Modeling Assumptions
- Lumped (0D) representation of plasma and wall inventories  
- Global flux balance governs plasma–wall interaction  
- Effective wall uptake and release rates (no spatial gradients)  
- Prescribed plasma exposure and external fueling  
- No explicit treatment of surface morphology, damage, or trap populations  

The model intentionally avoids spatial resolution and atomistic detail in order to remain transparent and computationally lightweight.

 Scope and Purpose
The purpose of this model is to capture dominant system-level hydrogen recycling behavior, such as wall inventory buildup and loss of density control, using minimal physics assumptions.

It is intended for conceptual analysis, sensitivity studies, and early-stage reasoning prior to the use of higher-fidelity plasma–wall interaction tools.

 Reference Case
The current parameter set is calibrated to WEST discharge conditions reported in:
Genesio et al., Nuclear Fusion 64 (2024).

This calibration serves as a reference case for evaluating model behavior rather than as a predictive validation.

 Status
Ongoing academic modeling effort.

