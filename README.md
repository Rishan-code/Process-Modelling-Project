# Phase-Field Simulation of Liposomal Nanocarrier Dissolution

## Project Overview
This project focuses on the computational modeling of liposomal nanocarrier dissolution and drug release under convective blood flow. The targeted delivery of therapeutics to the kidneys is hindered by rapid renal clearance and complex microvascular fluid dynamics. To address this, we employ a diffuse-interface Phase-Field approach (Convective Cahn-Hilliard equation) to naturally capture the moving boundary of the dissolving liposome, transitioning smoothly from the solid lipid core to the surrounding blood plasma. 

By balancing convective forces (Reynolds number, $Re$) and drug diffusion (Peclet number, $Pe$), this framework provides critical insights into shear-mediated drug release kinetics to inform the rational design of targeted renal nanomedicines.

## Progress So Far

### Phase 1: Problem Formulation and 1D Proof of Concept
- **Problem Statement & Governing Equations:** Defined the mathematical framework including the Convective Cahn-Hilliard equation for phase evolution and the mass balance equation for drug diffusion.
- **1D Implementation:** Developed a 1D Proof of Concept (`cahn_hilliard_1d.py`) to validate the theoretical hyperbolic tangent interface profile and establish the baseline for pure diffusion limiting cases.

### Phase 2: 2D Matrix Expansion and Convective Coupling
- **2D Mathematical Framework:** Expanded the model to 2D Cartesian coordinates to capture realistic flow phenomena around the spherical liposome.
- **Convective Flow Integration:** Integrated a 2D Poiseuille flow (parabolic velocity profile) to simulate laminar blood flow.
- **2D Implementation:** Developed a vectorized solver with heatmap output for the 2D model (`cahn_hilliard_2d.py`).

## Future Plans (Post-Midsem)
1. **Parametric Sweeps:** Systematically vary non-dimensional numbers, particularly $Re$ (characterizing flow velocity) and $Pe$ (characterizing advective to diffusive drug transport ratio), to evaluate their coupled effects on dissolution kinetics and release rates.
2. **Advanced Visualization:** Enhance transient 2D contour plots (heatmaps) of the phase variable and drug concentration, along with velocity streamlines, to visually elucidate the release mechanisms.
3. **Final Report Compilation:** Synthesize the methodology, parametric study results, and clinical implications into a comprehensive final report, formatted to graduate-level academic standards.

## Project Structure
- `phase1_deliverables.md`: Detailed documentation of Phase 1 including the problem statement, governing equations, and methodology.
- `phase2_deliverables.md`: Documentation of the 2D mathematical framework, initial conditions, and velocity field equations.
- `project_ledger.md`: Project activity and documentation ledger tracking tasks and status.
- `cahn_hilliard_1d.py`: Python script for the 1D phase-field model simulation.
- `cahn_hilliard_2d.py`: Python script for the 2D phase-field model simulation coupled with convective flow.
