 Reduced-Order Wall Recycling Model for Hydrogen in Fusion Reactor

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
 
Reference discharge context

A representative scenario is implemented for WEST discharge ("pulse") 57932.
Important: "calibrated to WEST pulse 57932" in this repository means:
- parameters were tuned here to reproduce reported qualitative features
  (timing/shape/order-of-magnitude behavior),
- not that every numeric parameter is directly copied from a single publication.

Parameter provenance rule
Every parameter must be classified as one of:
1) physical constant
2) chosen/assumed for demonstration
3) fitted in this repository
4) taken from literature (ONLY if explicitly cited and verifiable)
 Status
Ongoing academic modeling effort.

