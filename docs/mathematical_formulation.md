# Comprehensive Mathematical and Computational Formulation: Phase-Field Modeling of Liposomal Nanocarrier Dissolution and Drug Release

---

## 1. Physical System Description & Physiological Motivation
Targeted delivery of therapeutics to the kidneys (e.g., treating glomerulonephritis, polycystic kidney disease, or renal cell carcinoma) is constrained by two competing physiological barriers:
1. **Rapid Renal Clearance:** Unbound small-molecule drugs ($<6\text{ nm}$) are rapidly filtered out by the glomerular filtration barrier.
2. **Microvascular Hydrodynamic Shear:** In the renal microvasculature (arterioles and capillary networks), blood flow generates pronounced wall shear stresses ($\tau_w \approx 1 - 5\text{ Pa}$, shear rates $\dot{\gamma} \sim 10^2 - 10^3\text{ s}^{-1}$).

Liposomal nanocarriers ($R_0 \approx 50 - 100\text{ nm}$) protect encapsulated therapeutic payloads from premature glomerular elimination. However, upon entering the renal microcirculation, the liposomal lipid bilayer undergoes progressive degradation driven by thermodynamic undersaturation, enzymatic esterase cleavage, and hydrodynamic shear-mediated lipid shedding. 

Predicting the interfacial degradation and subsequent drug release requires a coupled continuum transport model that resolves:
* Moving boundaries and topological dissolution without explicit front tracking.
* Parabolic laminar blood flow (Poiseuille flow) in microchannels.
* Advection-diffusion of encapsulated solutes across a dynamic diffuse interface.

---

## 2. Phase-Field Thermodynamics & The Convective Cahn-Hilliard Equation

### 2.1 Order Parameter Definition
We define a continuous phase-field order parameter $\phi(\mathbf{x}, t) \in [0, 1]$:
$$\phi(\mathbf{x}, t) = \begin{cases} 0, & \text{inside the intact lipid bilayer/core (carrier phase)} \\ 1, & \text{in the surrounding blood plasma (continuous fluid phase)} \\ (0, 1), & \text{across the diffuse interfacial transition zone of characteristic thickness } \xi \end{cases}$$

### 2.2 Free Energy Functional
The total Ginzburg-Landau free energy functional $\mathcal{F}[\phi]$ of the heterogeneous system is formulated as:
$$\mathcal{F}[\phi] = \int_{\Omega} \left[ f_0(\phi) + \frac{\epsilon^2}{2} |\nabla \phi|^2 \right] d\Omega$$
where:
* $f_0(\phi)$ is the local bulk chemical free energy density, chosen as a canonical double-well potential with energy barrier height $W$:
  $$f_0(\phi) = W \phi^2 (1 - \phi)^2$$
  with two degenerate minima at $\phi = 0$ (lipid) and $\phi = 1$ (plasma).
* $\frac{\epsilon^2}{2} |\nabla \phi|^2$ represents the non-local gradient energy penalizing sharp compositional gradients, with gradient energy coefficient $\epsilon$.
* The equilibrium interfacial thickness $\xi$ and interfacial surface tension $\sigma$ are related by:
  $$\xi = \sqrt{2} \frac{\epsilon}{\sqrt{W}}, \qquad \sigma = \frac{\sqrt{2W}\epsilon}{6}$$

### 2.3 Variational Derivation of Chemical Potential
The chemical potential $\mu$ is defined as the functional (variational) derivative of $\mathcal{F}$ with respect to $\phi$:
$$\mu = \frac{\delta \mathcal{F}}{\delta \phi} = \frac{\partial f_0}{\partial \phi} - \nabla \cdot \left( \frac{\partial (\frac{\epsilon^2}{2} |\nabla \phi|^2)}{\partial (\nabla \phi)} \right) = \frac{\partial f_0}{\partial \phi} - \epsilon^2 \nabla^2 \phi$$
Evaluating the derivative of the double-well potential:
$$\frac{\partial f_0}{\partial \phi} = 2 W \phi (1 - \phi) (1 - 2\phi)$$
Thus:
$$\mu = 2 W \phi (1 - \phi) (1 - 2\phi) - \epsilon^2 \nabla^2 \phi$$

### 2.4 Degenerate Interfacial Mobility
To prevent unphysical diffusion or bulk Ostwald ripening within the pure phases and localize phase evolution strictly to the interfacial boundary layer, we define a degenerate mobility function:
$$M(\phi) = M_0 \phi^2 (1 - \phi)^2 + M_\infty$$
where $M_0$ is the characteristic mobility constant, and $M_\infty \ll M_0$ ($10^{-5} M_0$) is a numerical regularization term ensuring solver stability.

### 2.5 Thermodynamic Dissolution Driving Force & Shear Coupling
In physiological blood plasma, the lipid concentration is far below the critical micelle concentration (CMC), creating a thermodynamic driving force $\Delta \mu_{diss}$ promoting dissolution into the plasma phase ($\phi \to 1$). Furthermore, convective blood flow imposes hydrodynamic shear stress $\tau = \mu_{fluid} \dot{\gamma}$, which accelerates the stripping of surface lipid monomers. The generalized chemical potential becomes:
$$\mu = 2 W \phi (1 - \phi)(1 - 2\phi) - \Delta \mu_{diss} \cdot \Phi_{shear}(\dot{\gamma}) - \epsilon^2 \nabla^2 \phi$$
where $\Phi_{shear}(\dot{\gamma}) = 1 + \alpha_{shear} \frac{\dot{\gamma}}{\dot{\gamma}_0}$.

The convective Cahn-Hilliard equation governing the evolution of the liposome matrix is:
$$\frac{\partial \phi}{\partial t} + \mathbf{u} \cdot \nabla \phi = \nabla \cdot \left[ M(\phi) \nabla \mu \right]$$

---

## 3. Hydrodynamic Flow Field in Microchannels
Blood flow within microvascular capillaries and arterioles ($Re = \frac{\rho U_{max} H}{\mu_{fluid}} \ll 1$) is dominated by viscous forces. For a 2D channel of width $H$ ($y \in [0, H]$), the stationary velocity field $\mathbf{u} = (u_x, u_y)$ satisfies the Stokes equations and reduces to a fully developed Hagen-Poiseuille profile:
$$u_x(y) = 4 U_{max} \frac{y}{H} \left( 1 - \frac{y}{H} \right), \qquad u_y = 0$$
The corresponding local shear rate $\dot{\gamma}(y)$ is given by:
$$\dot{\gamma}(y) = \left| \frac{\partial u_x}{\partial y} \right| = \frac{4 U_{max}}{H} \left| 1 - \frac{2y}{H} \right|$$
This induces maximum shear stress $\tau_w = \mu_{fluid} \frac{4 U_{max}}{H}$ at the vessel boundaries ($y = 0, H$) and zero shear at the centerline ($y = H/2$).

---

## 4. Coupled Drug Mass Balance Equation

### 4.1 Transport & Phase-Dependent Diffusivity
The encapsulated therapeutic drug with local concentration $C(\mathbf{x}, t)$ is transported via advection and Fickian diffusion:
$$\frac{\partial C}{\partial t} + \mathbf{u} \cdot \nabla C = \nabla \cdot \left[ D(\phi) \nabla C \right] + S_{rel}(\phi, t)$$
The effective diffusion coefficient $D(\phi)$ varies over orders of magnitude between the dense hydrophobic lipid core ($D_{core}$) and the aqueous blood plasma ($D_{plasma} \gg D_{core}$):
$$D(\phi) = D_{core} + (D_{plasma} - D_{core}) h(\phi)$$
where $h(\phi)$ is a smooth monotonic interpolation function:
$$h(\phi) = \phi^2 (3 - 2\phi)$$
satisfying $h(0) = 0$, $h(1) = 1$, and $h'(0) = h'(1) = 0$.

### 4.2 Interfacial Release Source Term
As the lipid carrier erodes, encapsulated drug immobilized in the carrier is freed into the fluid phase. The volumetric release source rate is proportional to the local rate of matrix dissolution:
$$S_{rel}(\phi, \mathbf{x}, t) = k_{rel} (1 - \phi) C(\mathbf{x}, t) \cdot \mathbb{I}_{eroding}$$

---

## 5. Non-Dimensionalization & Dimensionless Groups

To generalize the simulation and identify the fundamental transport regimes, we introduce characteristic scaling scales:
* Characteristic length: Liposome initial radius $R_0$
* Characteristic velocity: Maximum centerline blood velocity $U_{max}$
* Characteristic convective time: $t_0 = \frac{R_0}{U_{max}}$
* Characteristic drug concentration: Initial loaded core concentration $C_0$
* Characteristic energy barrier: $W$

### 5.1 Dimensionless Variables
$$\mathbf{x}^* = \frac{\mathbf{x}}{R_0}, \quad \mathbf{u}^* = \frac{\mathbf{u}}{U_{max}}, \quad t^* = \frac{t U_{max}}{R_0}, \quad C^* = \frac{C}{C_0}, \quad \mu^* = \frac{\mu}{W}$$

### 5.2 Dimensionless Governing Equations
1. **Dimensionless Cahn-Hilliard Equation:**
   $$\frac{\partial \phi}{\partial t^*} + \mathbf{u}^* \cdot \nabla^* \phi = \frac{1}{Pe_\phi} \nabla^* \cdot \left[ M^*(\phi) \nabla^* \mu^* \right]$$
   where the non-dimensional chemical potential is:
   $$\mu^* = 2\phi(1-\phi)(1-2\phi) - \Delta \mu^* \cdot \left(1 + \beta \dot{\gamma}^*\right) - Cn^2 \nabla^{*2} \phi$$

2. **Dimensionless Drug Mass Balance Equation:**
   $$\frac{\partial C^*}{\partial t^*} + \mathbf{u}^* \cdot \nabla^* C^* = \frac{1}{Pe_c} \nabla^* \cdot \left[ D^*(\phi) \nabla^* C^* \right] + Da \cdot S^*(\phi, C^*)$$

### 5.3 Governing Dimensionless Groups
| Group | Formula | Physical Interpretation | Physiological Range |
| :--- | :--- | :--- | :--- |
| **Cahn Number ($Cn$)** | $\frac{\epsilon}{R_0}$ | Ratio of diffuse interface thickness to liposome radius | $0.02 - 0.08$ |
| **Péclet Number (Solute, $Pe_c$)** | $\frac{U_{max} R_0}{D_{plasma}}$ | Ratio of convective mass transport to molecular drug diffusion | $1 - 1000$ |
| **Phase Mobility Péclet ($Pe_\phi$)** | $\frac{U_{max} R_0^2}{M_0 W}$ | Ratio of advective timescale to interfacial relaxation timescale | $10 - 200$ |
| **Damköhler Number ($Da$)** | $\frac{k_{rel} R_0}{U_{max}}$ | Ratio of reaction/dissolution rate to convective flow rate | $0.01 - 5.0$ |
| **Diffusivity Contrast ($\delta_D$)** | $\frac{D_{core}}{D_{plasma}}$ | Ratio of lipid core diffusivity to bulk plasma diffusivity | $10^{-4} - 10^{-2}$ |
| **Sherwood Number ($Sh$)** | $\frac{k_m (2R_0)}{D_{plasma}}$ | Ratio of total mass transfer rate to pure diffusive conductance | Calculated from flux |

---

## 6. Numerical Discretization & Algorithmic Architecture

### 6.1 Spatial Discretization (Conservative Finite Difference)
* **Domain:** Rectangular 2D Cartesian grid $[0, L_x^*] \times [0, L_y^*]$ discretized with uniform spacing $\Delta x^*, \Delta y^*$.
* **Laplacian and Gradient:** Second-order central differences:
  $$\nabla^2 \phi_{i,j} = \frac{\phi_{i+1,j} - 2\phi_{i,j} + \phi_{i-1,j}}{(\Delta x)^2} + \frac{\phi_{i,j+1} - 2\phi_{i,j} + \phi_{i,j-1}}{(\Delta y)^2}$$
* **Conservative Flux Divergence:**
  $$\nabla \cdot \mathbf{J} = \frac{J_{x, i+1/2, j} - J_{x, i-1/2, j}}{\Delta x} + \frac{J_{y, i, j+1/2} - J_{y, i, j-1/2}}{\Delta y}$$
  where face mobilities and diffusivities are evaluated using arithmetic/harmonic mid-point averages to strictly ensure mass conservation.
* **Advection Scheme:** To eliminate spurious numerical oscillations without excessive numerical diffusion, advective terms $\mathbf{u} \cdot \nabla \phi$ and $\mathbf{u} \cdot \nabla C$ are discretized using a second-order flux-limited / upwind formulation:
  $$u_x \frac{\partial \phi}{\partial x} \approx u_x \left[ \frac{\phi_{i,j} - \phi_{i-1,j}}{\Delta x} \right] \quad (u_x \ge 0)$$

### 6.2 Boundary Conditions
* **Inlet ($x^* = 0$):**
  $$\phi(0, y^*, t^*) = 1, \quad C^*(0, y^*, t^*) = 0, \quad u_x^*(y^*) = 4 y^* (1 - y^* / L_y^*)$$
* **Outlet ($x^* = L_x^*$):**
  Convective outflow condition: $\frac{\partial \phi}{\partial x^*} = 0, \quad \frac{\partial C^*}{\partial x^*} = 0$
* **Walls ($y^* = 0, L_y^*$):**
  No-slip & zero-flux:
  $$u_x^* = 0, \quad \mathbf{n} \cdot \nabla \phi = 0, \quad \mathbf{n} \cdot \nabla \mu = 0, \quad \mathbf{n} \cdot \nabla C^* = 0$$

### 6.3 Mass Conservation Verification Formulation
The instantaneous total lipid mass $M_{lipid}(t)$ and retained drug mass $M_{drug}(t)$ are monitored via numerical quadrature:
$$M_{lipid}(t^*) = \iint_{\Omega} [1 - \phi(x^*, y^*, t^*)] \, dx^* dy^*$$
$$M_{drug}(t^*) = \iint_{\Omega} C^*(x^*, y^*, t^*) \, dx^* dy^*$$
The cumulative drug mass washed out across the outflow boundary is:
$$M_{out}(t^*) = \int_0^{t^*} \int_0^{L_y^*} u_x^*(y^*) C^*(L_x^*, y^*, \tau) \, dy^* d\tau$$
Total drug conservation requires:
$$M_{total}(t^*) = M_{drug}(t^*) + M_{out}(t^*) = M_{drug}(0) \quad (\pm \text{ error } < 0.5\%)$$

