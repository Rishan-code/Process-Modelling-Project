# Phase 2 Deliverables: 2D Phase-Field Simulation with Convective Flow

## 2D Mathematical Framework

### Convective Cahn-Hilliard Equation (2D Cartesian)
In two spatial dimensions $(x, y)$, the evolution of the phase parameter $\phi$ is governed by the 2D Convective Cahn-Hilliard equation:

$$ \frac{\partial \phi}{\partial t} + v_x \frac{\partial \phi}{\partial x} + v_y \frac{\partial \phi}{\partial y} = \nabla \cdot \left( M(\phi) \nabla \mu \right) $$

where the chemical potential $\mu$ in 2D is:

$$ \mu = \frac{d f(\phi)}{d \phi} - \epsilon^2 \left( \frac{\partial^2 \phi}{\partial x^2} + \frac{\partial^2 \phi}{\partial y^2} \right) $$

and $\mathbf{v} = (v_x, v_y)$ is the convective velocity vector. The mobility $M(\phi)$ remains a degenerate function $M_0 \phi^2(1-\phi)^2$ to restrict diffusion strictly to the interface.

### Drug Mass Balance Equation (2D Cartesian)
The transport of the drug concentration $C$ under convective flow and phase-dependent diffusion is given by:

$$ \frac{\partial C}{\partial t} + v_x \frac{\partial C}{\partial x} + v_y \frac{\partial C}{\partial y} = \nabla \cdot \left( D(\phi) \nabla C \right) + S(\phi, \mathbf{v}) $$

where $D(\phi)$ is the effective diffusivity. In 2D, the divergence operator expands as $\frac{\partial}{\partial x} \left( D(\phi) \frac{\partial C}{\partial x} \right) + \frac{\partial}{\partial y} \left( D(\phi) \frac{\partial C}{\partial y} \right)$.

### 2D Initial Conditions
To initialize a circular liposome of radius $R$ centered at $(x_c, y_c)$, we define the distance from the center $r = \sqrt{(x - x_c)^2 + (y - y_c)^2}$. The phase variable $\phi$ is initialized using a smoothed step function (hyperbolic tangent) to define the diffuse interface:

$$ \phi(x, y, t=0) = \frac{1}{2} \left[ 1 + \tanh\left( \frac{r - R}{\sqrt{2}\epsilon} \right) \right] $$

This formula guarantees $\phi \approx 0$ inside the solid lipid core ($r < R$), $\phi \approx 1$ in the bulk blood plasma ($r > R$), and a perfectly smooth continuous transition zone where $\phi = 0.5$ at the exact boundary ($r = R$).

### Velocity Field: 2D Poiseuille Flow
To simulate laminar blood flow in a vessel of width $H$ (spanning from $y = 0$ to $y = H$), we assume a fully developed 2D parabolic Poiseuille flow profile. The velocity vector $\mathbf{v} = (v_x, v_y)$ is strictly unidirectional in the longitudinal $x$-direction ($v_y = 0$):

$$ v_x(y) = 4 v_{max} \frac{y}{H} \left( 1 - \frac{y}{H} \right) $$

$$ v_y(y) = 0 $$

where $v_{max}$ is the maximum centerline velocity (at $y = H/2$). This stationary background flow field drives the convective transport within both the Cahn-Hilliard and mass balance equations.
