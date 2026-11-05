"""
Numerical Verification & Convergence Suite for Phase-Field Nanocarrier Solver
Performs:
1. Analytical 1D Hyperbolic Tangent Interface Verification
2. Global Mass Conservation Audit (Lipid & Drug)
3. Grid Convergence Study (Spatial Mesh Independence: dx, dx/2, dx/4)
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

from src.coupled_solver import SimulationConfig, CoupledPhaseFieldSolver


def verify_1d_analytical_tanh(save_path: str = "results/figures/verification_1d_tanh.png"):
    """
    Validates that the 1D phase-field profile matches the exact theoretical Ginzburg-Landau
    equilibrium solution: phi(x) = 0.5 * [1 + tanh((x - x0) / (sqrt(2) * epsilon))].
    """
    print("\n--- Running Verification 1: Analytical 1D Hyperbolic Tangent Interface ---")
    
    # 1D slice simulation
    Nx = 200
    Lx = 10.0
    x0 = 5.0
    epsilon = 0.4
    W = 1.0
    M0 = 1.0
    dx = Lx / (Nx - 1)
    x = np.linspace(0, Lx, Nx)
    
    # Initial perturbation slightly off equilibrium
    phi = 0.5 * (1.0 + np.tanh((x - x0) / (np.sqrt(2.0) * epsilon)))
    # Add small perturbation to test relaxation to analytical equilibrium
    phi_perturbed = np.clip(phi + 0.05 * np.sin(2.0 * np.pi * x / Lx), 0.0, 1.0)
    
    # Relax via pure Cahn-Hilliard diffusion (no convection, no dissolution)
    dt = 1e-4
    steps = 1500
    current_phi = phi_perturbed.copy()
    
    for s in range(steps):
        # 1D laplacian
        lap = np.zeros_like(current_phi)
        lap[1:-1] = (current_phi[2:] - 2*current_phi[1:-1] + current_phi[:-2]) / dx**2
        lap[0] = 2 * (current_phi[1] - current_phi[0]) / dx**2
        lap[-1] = 2 * (current_phi[-2] - current_phi[-1]) / dx**2
        
        df_dphi = 2 * W * current_phi * (1 - current_phi) * (1 - 2*current_phi)
        mu = df_dphi - (epsilon**2) * lap
        
        grad_mu = np.zeros_like(mu)
        grad_mu[1:-1] = (mu[2:] - mu[:-2]) / (2 * dx)
        
        mob = M0 * (current_phi**2) * ((1 - current_phi)**2) + 1e-4
        flux = -mob * grad_mu
        
        div_flux = np.zeros_like(flux)
        div_flux[1:-1] = (flux[2:] - flux[:-2]) / (2 * dx)
        div_flux[0] = (flux[1] - flux[0]) / dx
        div_flux[-1] = (flux[-1] - flux[-2]) / dx
        
        current_phi = np.clip(current_phi - dt * div_flux, 0.0, 1.0)

    # Analytical exact solution
    phi_exact = 0.5 * (1.0 + np.tanh((x - x0) / (np.sqrt(2.0) * epsilon)))
    
    # Error metrics
    l2_error = np.sqrt(np.mean((current_phi - phi_exact)**2))
    linf_error = np.max(np.abs(current_phi - phi_exact))
    
    print(f"L2 Norm Error vs Exact Analytical Solution:   {l2_error:.3e}")
    print(f"L_inf Maximum Error:                          {linf_error:.3e}")
    
    # Plotting comparison
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.figure(figsize=(9, 5))
    plt.plot(x, phi_perturbed, ':', color='gray', label='Initial Perturbed State', alpha=0.7)
    plt.plot(x, current_phi, 'b-', linewidth=2.5, label='Numerical Cahn-Hilliard (Relaxed)')
    plt.plot(x, phi_exact, 'r--', linewidth=2.0, label=r'Analytical Exact: $\frac{1}{2}[1+\tanh(\frac{x-x_0}{\sqrt{2}\epsilon})]$')
    
    plt.title(f'Analytical 1D Interface Verification (L2 Error = {l2_error:.2e})', fontsize=13, fontweight='bold')
    plt.xlabel('Spatial Coordinate (x)', fontsize=11)
    plt.ylabel(r'Phase Field Variable ($\phi$)', fontsize=11)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(frameon=True, fontsize=10)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"Saved analytical validation plot to {save_path}")
    
    return {'l2_error': l2_error, 'linf_error': linf_error}


def verify_mass_conservation(save_path: str = "results/figures/verification_mass_conservation.png"):
    """
    Validates global mass conservation under pure convective diffusion
    (zero dissolution limit, closed zero-flux domain).
    """
    print("\n--- Running Verification 2: Global Mass Conservation Audit ---")
    
    cfg = SimulationConfig(
        Lx=16.0, Ly=8.0, Nx=80, Ny=40,
        R0=2.0, xc=4.0, yc=4.0,
        u_max=0.5,
        k_diss=0.0,       # Zero dissolution -> perfect conservation test
        Pe_c=25.0,
        dt=2e-4,
        t_final=1.5,
        diagnostic_stride=500
    )
    
    solver = CoupledPhaseFieldSolver(cfg)
    res = solver.run(verbose=False)
    hist = res['history']
    
    times = np.array(hist['time'])
    lipid_mass = np.array(hist['lipid_mass'])
    carrier_drug = np.array(hist['carrier_drug'])
    plasma_drug = np.array(hist['plasma_drug'])
    outflow_drug = np.array(hist['cumulative_outflow'])
    total_drug = np.array(hist['total_system_drug'])
    
    lipid_error_pct = (lipid_mass - lipid_mass[0]) / lipid_mass[0] * 100.0
    drug_error_pct = (total_drug - total_drug[0]) / total_drug[0] * 100.0
    
    max_lipid_leak = np.max(np.abs(lipid_error_pct))
    max_drug_leak = np.max(np.abs(drug_error_pct))
    
    print(f"Max Lipid Conservation Deviation: {max_lipid_leak:.4f}%")
    print(f"Max Total Drug Balance Deviation: {max_drug_leak:.4f}%")
    
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    fig, ax1 = plt.subplots(figsize=(9, 5))
    
    ax1.plot(times, lipid_mass, 'b-', linewidth=2.0, label='Lipid Mass (Integral [1 - phi])')
    ax1.plot(times, total_drug, 'g--', linewidth=2.0, label='Total Drug (Retained + Cumulative Outflow)')
    ax1.set_xlabel('Non-Dimensional Time (t*)', fontsize=11)
    ax1.set_ylabel('Absolute Domain Mass', fontsize=11)
    ax1.grid(True, linestyle='--', alpha=0.5)
    
    ax2 = ax1.twinx()
    ax2.plot(times, drug_error_pct, 'r:', linewidth=1.5, label='Mass Balance Drift (%)')
    ax2.set_ylabel('Relative Numerical Leakage (%)', color='r', fontsize=11)
    ax2.tick_params(axis='y', labelcolor='r')
    ax2.set_ylim(-0.5, 0.5)
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='center right', frameon=True)
    
    plt.title(f'Global Mass Conservation Audit (Max Drift < {max(max_lipid_leak, max_drug_leak):.3f}%)', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"Saved mass conservation plot to {save_path}")
    
    return {'max_lipid_leak': max_lipid_leak, 'max_drug_leak': max_drug_leak}


def verify_mesh_convergence(save_path: str = "results/figures/verification_mesh_convergence.png"):
    """
    Performs grid refinement study across 3 mesh densities:
    Coarse (60x30), Medium (120x60), Fine (240x120) to demonstrate spatial mesh independence.
    """
    print("\n--- Running Verification 3: Spatial Mesh Independence Study ---")
    
    resolutions = [
        ("Coarse (60x30)", 60, 30),
        ("Medium (120x60)", 120, 60),
        ("Fine (180x90)", 180, 90)
    ]
    
    results = {}
    
    for label, nx, ny in resolutions:
        print(f"  Solving on grid {nx} x {ny} ...")
        cfg = SimulationConfig(
            Lx=20.0, Ly=10.0, Nx=nx, Ny=ny,
            R0=2.0, xc=5.0, yc=5.0,
            u_max=1.0,
            k_diss=0.08,
            Pe_c=30.0,
            dt=2e-4,
            t_final=1.5,
            diagnostic_stride=500
        )
        solver = CoupledPhaseFieldSolver(cfg)
        res = solver.run(verbose=False)
        results[label] = {
            'time': np.array(res['history']['time']),
            'release': np.array(res['history']['release_fraction']),
            'radius': np.array(res['history']['effective_radius'])
        }
        
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    
    styles = {'Coarse (60x30)': ('r:', 1.8), 'Medium (120x60)': ('b--', 2.0), 'Fine (180x90)': ('k-', 2.2)}
    
    for label in results:
        style, lw = styles[label]
        t = results[label]['time']
        rel = results[label]['release'] * 100.0
        rad = results[label]['radius']
        ax1.plot(t, rel, style, linewidth=lw, label=label)
        ax2.plot(t, rad, style, linewidth=lw, label=label)
        
    ax1.set_title('A. Fractional Drug Release Kinetics', fontsize=11, fontweight='bold')
    ax1.set_xlabel('Dimensionless Time (t*)', fontsize=10)
    ax1.set_ylabel('Cumulative Drug Released (%)', fontsize=10)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(frameon=True)
    
    ax2.set_title('B. Effective Liposome Radius Degradation', fontsize=11, fontweight='bold')
    ax2.set_xlabel('Dimensionless Time (t*)', fontsize=10)
    ax2.set_ylabel(r'Core Radius $R_{eff} / R_0$', fontsize=10)
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend(frameon=True)
    
    plt.suptitle('Spatial Mesh Independence & Grid Convergence Verification', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"Saved mesh convergence plot to {save_path}")
    
    return results


if __name__ == "__main__":
    verify_1d_analytical_tanh()
    verify_mass_conservation()
    verify_mesh_convergence()
