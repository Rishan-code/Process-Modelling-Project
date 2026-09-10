import numpy as np
import matplotlib.pyplot as plt

# Domain parameters
Lx, Ly = 20.0, 10.0      # Domain dimensions
Nx, Ny = 120, 60         # Number of grid points
dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)
x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(x, y, indexing='ij')

# Physical parameters
W = 1.0                  # Double-well energy barrier
epsilon = 0.4            # Gradient energy coefficient
M0 = 0.5                 # Baseline mobility
v_max = 3.0              # Maximum flow velocity

# Numerical parameters
dt = 2e-4                # Small time step for explicit stability
n_steps = 5000           # Total time steps
plot_freq = 1000         # Plot frequency

# Poiseuille velocity field (vx varies with y, vy = 0)
# Flow is in x-direction. y is from 0 to Ly.
vx = 4 * v_max * (Y / Ly) * (1 - Y / Ly)
vy = np.zeros_like(vx)

# Initial condition: Circular liposome
xc, yc = Lx / 4, Ly / 2  # Placed towards the left to allow space for flow
R = 2.5                  # Liposome radius
r = np.sqrt((X - xc)**2 + (Y - yc)**2)
phi = 0.5 * (1 + np.tanh((r - R) / (np.sqrt(2) * epsilon)))

# Fast Finite Difference Operators (Vectorized)
def laplacian_2d(p, dx, dy):
    lap = np.zeros_like(p)
    lap[1:-1, 1:-1] = (
        (p[2:, 1:-1] - 2*p[1:-1, 1:-1] + p[:-2, 1:-1]) / dx**2 +
        (p[1:-1, 2:] - 2*p[1:-1, 1:-1] + p[1:-1, :-2]) / dy**2
    )
    # Neumann BCs (zero derivative) at boundaries
    lap[0, :] = lap[1, :]
    lap[-1, :] = lap[-2, :]
    lap[:, 0] = lap[:, 1]
    lap[:, -1] = lap[:, -2]
    return lap

def grad_2d(p, dx, dy):
    gx = np.zeros_like(p)
    gy = np.zeros_like(p)
    gx[1:-1, :] = (p[2:, :] - p[:-2, :]) / (2*dx)
    gy[:, 1:-1] = (p[:, 2:] - p[:, :-2]) / (2*dy)
    return gx, gy

def div_2d(fx, fy, dx, dy):
    div = np.zeros_like(fx)
    div[1:-1, :] += (fx[2:, :] - fx[:-2, :]) / (2*dx)
    div[:, 1:-1] += (fy[:, 2:] - fy[:, :-2]) / (2*dy)
    return div

# Upwind scheme for advection to ensure stability in convective regime
def advection_upwind(p, vx, dx):
    adv_x = np.zeros_like(p)
    # Since vx >= 0 everywhere, use backward difference
    adv_x[1:, :] = vx[1:, :] * (p[1:, :] - p[:-1, :]) / dx
    adv_x[0, :] = vx[0, :] * (p[1, :] - p[0, :]) / dx
    return adv_x

if __name__ == "__main__":
    print("Starting 2D Convective Cahn-Hilliard simulation...")

    plt.figure(figsize=(15, 8))
    plot_idx = 1

    for step in range(n_steps + 1):
        # 1. Chemical potential
        df_dphi = 2 * W * phi * (1 - phi) * (1 - 2 * phi)
        mu = df_dphi - (epsilon**2) * laplacian_2d(phi, dx, dy)
        
        # 2. Degenerate Mobility
        mob = M0 * (phi**2) * ((1 - phi)**2) + 1e-4
        
        # 3. Flux
        g_mu_x, g_mu_y = grad_2d(mu, dx, dy)
        jx = -mob * g_mu_x
        jy = -mob * g_mu_y
        
        # 4. Phase-field evolution
        div_j = div_2d(jx, jy, dx, dy)
        adv = advection_upwind(phi, vx, dx)
        
        phi = phi + dt * (-div_j - adv)
        
        # Store heatmaps
        if step % plot_freq == 0:
            print(f"Step {step}/{n_steps} completed.")
            plt.subplot(2, 3, plot_idx)
            plt.pcolormesh(X, Y, phi, cmap='coolwarm', shading='gouraud', vmin=0, vmax=1)
            plt.colorbar(label=r'$\phi$')
            plt.title(f'Time: {step * dt:.2f}')
            plt.xlabel('x')
            plt.ylabel('y')
            plot_idx += 1

    plt.tight_layout()
    plt.savefig('cahn_hilliard_2d_result.png')
    print("Simulation complete. 2D spatiotemporal heatmap saved as 'cahn_hilliard_2d_result.png'")
