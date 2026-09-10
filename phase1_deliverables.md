# Phase 1 Deliverables: Phase-Field Simulation of Liposomal Nanocarrier Dissolution

## Problem Statement
The targeted delivery of therapeutics to the kidneys is significantly hindered by rapid renal clearance and complex microvascular fluid dynamics. Liposomal nanocarriers present a promising vehicle to protect hydrophobic drug payloads and facilitate controlled release. However, predicting their dissolution kinetics under physiological convective blood flow remains a formidable engineering challenge. Traditional sharp-interface tracking methods (e.g., Level Set or Volume of Fluid) are computationally expensive and struggle with topological changes during liposome degradation. To overcome this, we employ a diffuse-interface Phase-Field approach. This method naturally captures the moving boundary of the dissolving liposome—transitioning smoothly from the solid lipid core ($\phi = 0$) to the surrounding blood plasma ($\phi = 1$)—without explicit boundary tracking. By balancing convective forces (Reynolds number, $Re$) and drug diffusion (Peclet number, $Pe$), this framework will provide critical insights into shear-mediated drug release kinetics, ultimately informing the rational design of targeted renal nanomedicines.

## Schematic of the Model
**Text Prompt for Diagram Generation:**
"Create a technical scientific diagram showing a spherical liposomal nanocarrier dissolving in a bloodstream. On the left side of the liposome, depict convective blood flow using horizontal arrows labeled with Reynolds number '$Re$'. Surrounding the liposome, depict smaller, radiating arrows representing drug diffusion, labeled with Peclet number '$Pe$'. The liposome itself should have a distinct, solid core (labeled 'Solid Lipid Core, $\phi = 0$') that transitions into a blurry, gradient-filled outer corona (labeled 'Diffuse Interface'). The surrounding fluid should be labeled 'Blood Plasma, $\phi = 1$'. Below this physical representation, include a 2D line graph plotting the phase parameter '$\phi$' on the y-axis (from 0 to 1) against the spatial coordinate '$x$' on the x-axis, showing a smooth sigmoidal curve transitioning from $\phi=0$ inside the core to $\phi=1$ in the plasma."

## Governing Equations

The evolution of the liposome interface is governed by the Convective Cahn-Hilliard equation:

$$ \frac{\partial \phi}{\partial t} + \mathbf{u} \cdot \nabla \phi = \nabla \cdot \left( M(\phi) \nabla \mu \right) $$

where:
*   $\phi$ is the phase-field variable ($\phi \in [0, 1]$).
*   $\mathbf{u}$ is the convective velocity vector of the blood flow.
*   $M(\phi)$ is the phase-dependent mobility function.
*   $\mu$ is the chemical potential, defined as the variational derivative of the free energy functional $F$:

$$ \mu = \frac{\delta F}{\delta \phi} = \frac{d f(\phi)}{d \phi} - \epsilon^2 \nabla^2 \phi $$

Here, $f(\phi) = W \phi^2 (1-\phi)^2$ is the double-well bulk energy density (with barrier height $W$), and $\epsilon$ is the gradient energy coefficient controlling the interface thickness.

To ensure that dissolution dynamics (mobility) are strictly localized to the diffuse interface and vanish in the bulk phases, we define a degenerate mobility function:

$$ M(\phi) = M_0 \phi^2 (1-\phi)^2 $$

where $M_0$ is a constant baseline mobility.

The transient mass diffusion of the drug (concentration $C$) is coupled to the phase field:

$$ \frac{\partial C}{\partial t} + \mathbf{u} \cdot \nabla C = \nabla \cdot \left( D(\phi) \nabla C \right) + S(\phi, \mathbf{u}) $$

where $D(\phi)$ is the effective diffusion coefficient (which can be interpolated between the solid core and plasma, e.g., $D(\phi) = D_s(1-\phi) + D_p\phi$), and $S$ represents a source term coupled to the dissolution of the lipid matrix. 

## Solution Methodology
The system of coupled, non-linear partial differential equations will be solved using a continuum-level Eulerian approach. 
*   **Spatial Discretization:** The domain will be discretized using an Eulerian grid. Standard finite difference methods (FDM) will be employed to approximate the spatial derivatives (gradient, divergence, and Laplacian operators). Upwind schemes will be utilized for the convective terms to maintain numerical stability in advection-dominated regimes (high $Re$).
*   **Time Integration:** The Cahn-Hilliard equation involves fourth-order spatial derivatives (through the chemical potential), introducing severe stiffness to the system. Therefore, temporal integration will be handled using implicit, stable stiff ODE solvers. In Python, we will leverage `scipy.integrate.solve_ivp` utilizing the Backward Differentiation Formula (BDF) method, or alternatively, employ the `FiPy` library which natively handles implicit solutions for phase-field problems.

## Validation Strategy
The numerical model will be rigorously validated through several fundamental checks:
1.  **Limiting Case - Pure Diffusion:** By setting the velocity field $\mathbf{u} = 0$ ($Re \to 0$), the model's drug release profile should collapse to the analytical solutions for pure Fickian diffusion from a sphere.
2.  **Mass Conservation:** The total mass of the system (integral of $\phi$ and $C$ over the entire domain) must be strictly conserved across the phase boundary, assuming zero-flux boundary conditions on the outer domain walls.
3.  **Equilibrium Interface Profile:** In the absence of flow and dissolution driving forces, the 1D phase profile across the interface must precisely match the theoretical hyperbolic tangent solution: $\phi(x) = \frac{1}{2} \left[ 1 + \tanh\left(\frac{x}{\sqrt{2}\epsilon}\right) \right]$.

## Plan for Post-Midsem
1.  **Transition to 2D Spatial Matrices:** Expand the 1D model to a 2D Cartesian or axisymmetric grid to capture realistic flow phenomena around the spherical liposome (e.g., wake formation, shear stress gradients).
2.  **Parametric Sweep of Non-Dimensional Numbers:** Systematically vary $Re$ (characterizing flow velocity) and $Pe$ (characterizing the ratio of advective to diffusive drug transport) to evaluate their coupled effects on dissolution kinetics and release rates.
3.  **Visualization and Heatmap Generation:** Generate transient 2D contour plots (heatmaps) of the phase variable $\phi$ and drug concentration $C$, along with velocity streamlines, to visually elucidate the release mechanisms.
4.  **Final Report Compilation:** Synthesize the methodology, parametric study results, and clinical implications into a comprehensive final report, formatted to graduate-level academic standards.
