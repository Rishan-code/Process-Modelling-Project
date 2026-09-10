import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parameters
L = 10.0         # Domain length
N = 100          # Number of grid points
dx = L / (N - 1) # Grid spacing
x = np.linspace(0, L, N)

W = 1.0          # Double-well energy barrier
epsilon = 0.5    # Gradient energy coefficient
M0 = 1.0         # Mobility constant

# Degenerate mobility function
def M(phi):
    return M0 * (phi**2) * ((1 - phi)**2)

# Derivative of the double well potential
def df_dphi(phi):
    return 2 * W * phi * (1 - phi) * (1 - 2 * phi)

# Spatial derivatives (using central differences with zero-flux Neumann boundary conditions)
def laplacian(phi, dx):
    lap = np.zeros_like(phi)
    lap[1:-1] = (phi[2:] - 2*phi[1:-1] + phi[:-2]) / dx**2
    lap[0] = 2 * (phi[1] - phi[0]) / dx**2
    lap[-1] = 2 * (phi[-2] - phi[-1]) / dx**2
    return lap

def gradient(phi, dx):
    grad = np.zeros_like(phi)
    grad[1:-1] = (phi[2:] - phi[:-2]) / (2 * dx)
    grad[0] = 0
    grad[-1] = 0
    return grad

def divergence(flux, dx):
    div = np.zeros_like(flux)
    div[1:-1] = (flux[2:] - flux[:-2]) / (2 * dx)
    # Simple one-sided differences at boundaries for divergence
    div[0] = (flux[1] - flux[0]) / dx
    div[-1] = (flux[-1] - flux[-2]) / dx
    return div

# Cahn-Hilliard equation (RHS for ODE solver)
def cahn_hilliard_1d(t, phi):
    # Chemical potential: mu = df/dphi - epsilon^2 * Laplacian(phi)
    mu = df_dphi(phi) - (epsilon**2) * laplacian(phi, dx)
    
    # Mobility
    mob = M(phi)
    # Add a small offset to mobility to prevent it from being exactly zero initially
    # which can stall the solver if the interface is perfectly sharp.
    mob = mob + 1e-4 
    
    # Flux: J = -M(phi) * grad(mu)
    grad_mu = gradient(mu, dx)
    flux = -mob * grad_mu
    
    # dphi/dt = -div(Flux)
    dphi_dt = -divergence(flux, dx)
    return dphi_dt

# Initial condition: Step function (sharp interface)
phi0 = np.zeros(N)
phi0[x > L/2] = 1.0

# Time span for simulation
t_span = (0, 2.0)
t_eval = np.linspace(0, 2.0, 5)

if __name__ == "__main__":
    print("Solving 1D Cahn-Hilliard equation. This may take a moment...")
    # Solve the PDE using the method of lines
    sol = solve_ivp(cahn_hilliard_1d, t_span, phi0, method='BDF', t_eval=t_eval, vectorized=False)

    # Plotting the results
    plt.figure(figsize=(10, 6))
    for i, t in enumerate(sol.t):
        plt.plot(x, sol.y[:, i], label=f't = {t:.2f}', marker='.' if i==0 else '')

    plt.title('1D Cahn-Hilliard: Interface Smoothing (Proof of Concept)')
    plt.xlabel('Spatial Coordinate (x)')
    plt.ylabel(r'Phase Parameter ($\phi$)')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('cahn_hilliard_1d_result.png')
    print("Simulation complete. Plot saved as 'cahn_hilliard_1d_result.png'")
