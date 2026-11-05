# Phase-Field Simulation of Liposomal Nanocarrier Dissolution & Drug Release

[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![NumPy](https://img.shields.io/badge/NumPy-Vectorized%20FDM-brightgreen.svg)](https://numpy.org/)
[![SciPy](https://img.shields.io/badge/SciPy-PDE%20%26%20Optimization-orange.svg)](https://scipy.org/)
[![Status](https://img.shields.io/badge/Status-Complete%20Capstone-success.svg)](#)

A comprehensive continuum-level computational transport framework modeling the dissolution kinetics and drug release of liposomal nanocarriers under convective blood flow in renal microvasculature.

---

## 🌟 Executive Highlights

* **Diffuse-Interface Phase-Field Dynamics:** Formulated and solved the 2D Convective Cahn-Hilliard equation with degenerate interfacial mobility, tracking the moving boundary of the dissolving nanocarrier without sharp-interface reconstruction singularities.
* **Multiphysics Coupling:** Simultaneously resolves:
  1. Fourth-order Cahn-Hilliard phase evolution ($\phi$).
  2. Hagen-Poiseuille parabolic laminar blood flow ($u_x(y)$) and wall shear rate ($\dot{\gamma}(y)$).
  3. Interfacial shear-enhanced matrix degradation ($r_{diss} \propto \phi(1-\phi)(1 + \beta_{shear}\dot{\gamma}^*)$).
  4. Convective-diffusion transport of encapsulated therapeutics ($C$) with phase-dependent tensor diffusivity ($D_{plasma} / D_{core} = 10^2$).
* **Numerical Verification & Stability:**
  * Exact match against analytical 1D hyperbolic tangent equilibrium interface ($L_2\text{ error} = 3.52 \times 10^{-2}$).
  * Rigorous mass conservation audit across closed domains (**<0.021% numerical leakage** over 15,000 steps).
  * Spatial grid independence demonstrated across $60\times30$, $120\times60$, and $180\times90$ meshes (<1.2% variation).
* **Levich Boundary-Layer Scaling:** Demonstrated that the convective mass transfer Sherwood number scales as **$Sh \propto Pe^{0.341}$**, matching theoretical Levich creeping shear boundary-layer theory ($Sh \sim Pe^{1/3}$) within 2.4% error.
* **Pharmacokinetic Model Fitting:** Release kinetics confirm anomalous non-Fickian transport (**Korsmeyer-Peppas exponent $n = 0.624$, $R^2 = 0.991$**), co-governed by Fickian diffusion and shear-mediated matrix erosion.

---

## 📁 Repository Structure

```
Process Modelling Project/
├── src/
│   ├── __init__.py
│   ├── coupled_solver.py        # Core 2D coupled solver (Cahn-Hilliard + Flow + Solute + Shear)
│   ├── verification.py          # 1D analytical tanh, mass conservation & mesh convergence suite
│   ├── parametric_sweeps.py     # Automated sweeps over Pe, beta_shear, and Cahn number (Cn)
│   └── kinetics_analysis.py     # Pharmacokinetic model fitting & Sherwood (Sh vs Pe) scaling
├── results/
│   ├── figures/                 # Publication-quality figures & simulation snapshots
│   │   ├── coupled_simulation_summary.png
│   │   ├── verification_1d_tanh.png
│   │   ├── verification_mass_conservation.png
│   │   ├── verification_mesh_convergence.png
│   │   ├── parametric_peclet_release.png
│   │   ├── parametric_peclet_plumes.png
│   │   ├── parametric_shear_coupling.png
│   │   ├── parametric_cahn_number.png
│   │   ├── pharmacokinetics_models.png
│   │   └── sherwood_vs_peclet.png
│   └── data/                    # Exported metrics & kinetic summary tables
├── docs/
│   ├── FINAL_REPORT.md          # Comprehensive graduate-level technical thesis (~35 pages caliber)
│   ├── mathematical_formulation.md # Complete thermodynamic derivations & non-dimensionalization
│   ├── phase1_deliverables.md   # Midterm Phase 1 problem formulation & 1D PoC
│   └── phase2_deliverables.md   # Midterm Phase 2 2D framework & Poiseuille coupling
├── project_ledger.md            # Complete 14-week semester activity and milestones ledger
└── README.md                    # Project overview & documentation
```

---

## 🚀 Quickstart & Execution

All simulation scripts are modular and can be executed directly from the root project directory:

### 1. Run the Production Coupled 2D Simulation
Solves the coupled system and generates the 4-panel multi-field summary plot:
```bash
python -m src.coupled_solver
```

### 2. Execute the Numerical Verification Suite
Validates the analytical hyperbolic-tangent profile, audits global mass conservation, and verifies mesh independence:
```bash
python -m src.verification
```

### 3. Run the Multi-Dimensional Parametric Sweep
Simulates across Péclet numbers ($Pe \in [5, 400]$), shear couplings, and Cahn numbers:
```bash
python -m src.parametric_sweeps
```

### 4. Pharmacokinetic Modeling & Sherwood Scaling
Fits Higuchi, Korsmeyer-Peppas, and Hixson-Crowell models and extracts $Sh(Pe)$ boundary layer scaling:
```bash
python -m src.kinetics_analysis
```

---

## 🔬 Core Results & Visualizations

### 1. Coupled Multiphysics Field Evolution
![Coupled Simulation Summary](results/figures/coupled_simulation_summary.png)
*(A) Phase field $\phi$ moving interface; (B) Hagen-Poiseuille streamlines; (C) Downstream drug concentration plume; (D) Core radius degradation and fractional release kinetics.*

### 2. Numerical Verification Suite
| 1D Analytical Tanh Profile | Global Mass Conservation Audit | Grid Independence Verification |
| :---: | :---: | :---: |
| ![1D Tanh](results/figures/verification_1d_tanh.png) | ![Mass Conservation](results/figures/verification_mass_conservation.png) | ![Mesh Convergence](results/figures/verification_mesh_convergence.png) |
| $L_2\text{ error} = 3.52 \times 10^{-2}$ | Numerical leakage < 0.021% | Coarse, Medium & Fine grids collapse |

### 3. Transport Scaling & Pharmacokinetic Kinetics
| Sherwood vs Péclet Scaling ($Sh \propto Pe^{1/3}$) | Empirical Pharmaceutical Model Fitting |
| :---: | :---: |
| ![Sherwood Scaling](results/figures/sherwood_vs_peclet.png) | ![Pharmacokinetics](results/figures/pharmacokinetics_models.png) |
| Matches Levich boundary-layer scaling within 2.4% | Korsmeyer-Peppas $R^2 = 0.991$, $n = 0.624$ |

---

## 📖 Complete Documentation
* For complete thermodynamic derivations and mathematical formulation, see [mathematical_formulation.md](docs/mathematical_formulation.md).
* For the comprehensive thesis-style capstone report, see [FINAL_REPORT.md](docs/FINAL_REPORT.md).
* For the full 14-week semester progression, see [project_ledger.md](project_ledger.md).
