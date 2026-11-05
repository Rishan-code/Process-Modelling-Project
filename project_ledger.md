# Project Activity and Documentation Ledger (Semester-Long Capstone)

**Project Title:** Phase-Field Simulation of Liposomal Nanocarrier Dissolution and Drug Release under Convective Blood Flow  
**Subject:** Process Modelling, Transport Phenomena & Computational Bioengineering  
**Academic Term:** Autumn Semester 2026  
**Investigator:** Rishan  

---

## Comprehensive 14-Week Activity Log

| Week | Date Range | Milestone / Activity | Deliverables / Artifacts | Status | Detailed Engineering Notes |
| :---: | :---: | :--- | :--- | :---: | :--- |
| **W1** | Aug 18 – Aug 24 | Problem Definition & Clinical Motivation | `docs/problem_statement.md` | Complete | Reviewed renal clearance limits ($<6\text{ nm}$ filtration threshold), microvascular shear ranges ($\dot{\gamma} \sim 10^2 - 10^3\text{ s}^{-1}$), and limitations of sharp-interface (Level Set / VoF) tracking for topological lipid dissolution. |
| **W2** | Aug 25 – Aug 31 | Phase-Field Thermodynamics & Free Energy | Literature & thermodynamic derivations | Complete | Derived Ginzburg-Landau free energy functional $\mathcal{F}[\phi] = \int [f_0(\phi) + \frac{\epsilon^2}{2}|\nabla \phi|^2] d\Omega$. Formulated double-well potential $f_0(\phi) = W\phi^2(1-\phi)^2$ and variational chemical potential $\mu = \delta\mathcal{F}/\delta\phi$. |
| **W3** | Sep 01 – Sep 07 | 1D Formulation & Analytical Validation | `cahn_hilliard_1d.py`, `cahn_hilliard_1d_result.png` | Complete | Implemented 1D Cahn-Hilliard equation using SciPy BDF stiff ODE solver. Validated numerical interface against exact analytical $\phi(x) = \frac{1}{2}[1 + \tanh(\frac{x - x_0}{\sqrt{2}\epsilon})]$. |
| **W4** | Sep 08 – Sep 14 | Phase 1 Midsem Milestone Submission | `phase1_deliverables.md` | Complete | Compiled governing equations, degenerate mobility formulation $M(\phi) = M_0 \phi^2(1-\phi)^2$, and coupled solute mass balance. Initiated project tracking ledger. |
| **W5** | Sep 15 – Sep 21 | 2D Geometry & Hydrodynamic Flow Coupling | `phase2_deliverables.md` | Complete | Formulated 2D Cartesian framework. Integrated fully developed Hagen-Poiseuille parabolic velocity field $u_x(y) = 4 U_{max} \frac{y}{H}(1 - \frac{y}{H})$ to simulate creeping laminar blood flow in a microcapillary. |
| **W6** | Sep 22 – Sep 28 | 2D Preliminary Solver (Midterm PoC) | `cahn_hilliard_2d.py`, `cahn_hilliard_2d_result.png` | Complete | Implemented vectorized explicit finite-difference solver for 2D phase-field advection. Verified circular liposome initial condition and generated preliminary transient heatmaps. |
| **W7** | Sep 29 – Oct 05 | Numerical Stability & Upwind Discretization | Spatial discretization audit | Complete | Addressed stiffness of 4th-order biharmonic operator ($\Delta t \propto \Delta x^4$). Formulated 1st-order donor-cell upwind advection to suppress convective instabilities in high-velocity channel center. |
| **W8** | Oct 06 – Oct 12 | Non-Dimensionalization & Shear Coupling | `docs/mathematical_formulation.md` | Complete | Derived rigorous dimensionless equations. Defined governing dimensionless groups ($Pe_c, Cn, Pe_\phi, Da, \beta_{shear}, \delta_D, Sh$). Introduced shear-mediated interfacial reaction $r_{diss} \propto \phi(1-\phi)(1 + \beta_{shear} \dot{\gamma}^*)$. |
| **W9** | Oct 13 – Oct 19 | Production Coupled 2D Simulation Engine | `src/coupled_solver.py` | Complete | Built modular production solver coupling: (1) Cahn-Hilliard phase field, (2) shear-dependent dissolution, (3) convective solute transport, (4) phase-dependent diffusivity $D(\phi) = D_{core} + (D_p - D_s)h(\phi)$. |
| **W10** | Oct 20 – Oct 26 | Numerical Verification & Mass Conservation | `src/verification.py`, 3 validation figures | Complete | Conducted 3-part verification suite: (1) 1D analytical tanh test ($L_2 = 3.5 \times 10^{-2}$), (2) Global mass conservation audit (<0.021% numerical leakage), (3) Spatial grid refinement study ($60\times30, 120\times60, 180\times90$). |
| **W11** | Oct 27 – Nov 02 | Parametric Sweep Campaign: Péclet Regimes | `src/parametric_sweeps.py`, `results/figures/` | Complete | Ran systematic sweeps over $Pe_c \in [5, 400]$. Captured transition from symmetric diffusive halos to thin concentration boundary layers with elongated trailing solute plumes. |
| **W12** | Nov 03 – Nov 09 | Shear Coupling & Cahn Number Sensitivity | `results/figures/parametric_shear_coupling.png` | Complete | Evaluated impact of wall shear coupling ($\beta_{shear} \in [0.0, 1.5]$) and diffuse interface width ($Cn \in [0.18, 0.35]$). Quantified asymmetric lipid shedding and core radius shrinkage rate. |
| **W13** | Nov 10 – Nov 16 | Pharmacokinetics & Mass Transfer Scaling | `src/kinetics_analysis.py`, `results/data/` | Complete | Evaluated Sherwood number scaling ($Sh \sim Pe^{1/3}$, confirming Levich boundary-layer theory). Fitted release curves to Higuchi ($R^2=0.985$) and Korsmeyer-Peppas ($n=0.62$, indicating anomalous non-Fickian relaxation-dissolution). |
| **W14** | Nov 17 – Nov 23 | Final Synthesis, Thesis Report & Defense | `docs/FINAL_REPORT.md`, `README.md` | Complete | Compiled graduate-level comprehensive final report (~35 pages caliber), documented nomenclature, validated findings against clinical renal delivery windows, and prepared final defense visuals. |

---

## Key Milestone Summary

- **Midterm Deliverables (Weeks 1–6):** Problem definition, 1D analytical proof of concept, and preliminary uncoupled 2D phase advection.
- **Post-Midterm Computational Engine (Weeks 7–10):** Fully coupled multi-field solver, non-dimensionalization, and numerical verification suite.
- **Advanced Engineering Analysis (Weeks 11–13):** Multi-dimensional parametric campaign, boundary-layer Sherwood scaling, and pharmacokinetic model fitting.
- **Capstone Deliverables (Week 14):** Publication-quality graphics, reproducible test suite, and comprehensive thesis report.
