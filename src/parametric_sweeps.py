"""
Automated Parametric Sweep Suite for Coupled Phase-Field Nanocarrier Dissolution
Systematically varies:
1. Péclet Number (Pe_c in [5, 20, 50, 150, 400]) -> Advection vs Diffusion
2. Shear Coupling / Dissolution Rate (k_diss and beta_shear) -> Hydrodynamic erosion
3. Nanocarrier Size / Cahn Number (Cn in [0.15, 0.25, 0.40]) -> Geometry scaling
Outputs high-resolution multi-panel comparison figures and data tables.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import time

from src.coupled_solver import SimulationConfig, CoupledPhaseFieldSolver


def run_peclet_sweep(pe_values=[5.0, 20.0, 50.0, 150.0, 400.0], save_dir="results"):
    """Evaluates the effect of Péclet number on drug plume dynamics and release rate."""
    print("\n========================================================")
    print("  PARAMETRIC SWEEP 1: PÉCLET NUMBER (Pe_c) VARIATION")
    print("========================================================")
    
    os.makedirs(f"{save_dir}/figures", exist_ok=True)
    os.makedirs(f"{save_dir}/data", exist_ok=True)
    
    results = {}
    
    for pe in pe_values:
        print(f"--> Simulating Pe_c = {pe:5.1f} ...")
        cfg = SimulationConfig(
            Lx=20.0, Ly=10.0, Nx=100, Ny=50,
            R0=2.2, xc=5.0, yc=5.0,
            u_max=1.5,
            k_diss=0.08,
            Pe_c=pe,
            dt=2e-4,
            t_final=2.5,
            diagnostic_stride=500
        )
        solver = CoupledPhaseFieldSolver(cfg)
        res = solver.run(verbose=False)
        
        results[pe] = {
            'time': np.array(res['history']['time']),
            'release': np.array(res['history']['release_fraction']),
            'lipid_mass': np.array(res['history']['lipid_mass']),
            'effective_radius': np.array(res['history']['effective_radius']),
            'final_C': res['C'].copy(),
            'X': res['X'].copy(),
            'Y': res['Y'].copy()
        }
    
    # 1. Multi-panel Release Curves Comparison
    plt.figure(figsize=(10, 6))
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(pe_values)))
    
    for i, pe in enumerate(pe_values):
        t = results[pe]['time']
        rel = results[pe]['release'] * 100.0
        plt.plot(t, rel, label=f'Pe = {pe:.0f}', color=colors[i], linewidth=2.2)
        
    plt.title('Effect of Péclet Number on Transient Drug Release Kinetics', fontsize=13, fontweight='bold')
    plt.xlabel('Dimensionless Time (t*)', fontsize=11)
    plt.ylabel('Fractional Drug Released (%)', fontsize=11)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(frameon=True, fontsize=10, loc='lower right')
    plt.tight_layout()
    plt.savefig(f"{save_dir}/figures/parametric_peclet_release.png", dpi=300)
    plt.close()
    
    # 2. 2D Concentration Plume Heatmaps across Pe
    fig, axes = plt.subplots(len(pe_values), 1, figsize=(11, 2.3 * len(pe_values)), sharex=True)
    if len(pe_values) == 1:
        axes = [axes]
        
    for i, pe in enumerate(pe_values):
        ax = axes[i]
        c = ax.pcolormesh(results[pe]['X'], results[pe]['Y'], results[pe]['final_C'], 
                          cmap='magma', shading='gouraud', vmin=0, vmax=1)
        ax.set_ylabel(f'Pe={pe:.0f}\ny / R0', fontsize=10, fontweight='bold')
        if i == len(pe_values) - 1:
            ax.set_xlabel('Channel Length (x / R0)', fontsize=11)
        if i == 0:
            ax.set_title(f'Spatial Drug Plume Profiles at t* = 2.5 across Péclet Regimes', fontsize=12, fontweight='bold')
            
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(c, cax=cbar_ax, label=r'Concentration $C / C_0$')
    plt.savefig(f"{save_dir}/figures/parametric_peclet_plumes.png", dpi=300)
    plt.close()
    
    print(f"Completed Péclet sweep. Figures saved to {save_dir}/figures/")
    return results


def run_shear_sweep(beta_values=[0.0, 0.4, 0.8, 1.5], save_dir="results"):
    """Evaluates the effect of hydrodynamic shear-induced erosion on matrix degradation."""
    print("\n========================================================")
    print("  PARAMETRIC SWEEP 2: SHEAR COUPLING (beta_shear) VARIATION")
    print("========================================================")
    
    os.makedirs(f"{save_dir}/figures", exist_ok=True)
    results = {}
    
    for beta in beta_values:
        print(f"--> Simulating beta_shear = {beta:4.2f} ...")
        cfg = SimulationConfig(
            Lx=20.0, Ly=10.0, Nx=100, Ny=50,
            R0=2.2, xc=5.0, yc=5.0,
            u_max=1.5,
            k_diss=0.08,
            beta_shear=beta,
            Pe_c=30.0,
            dt=2e-4,
            t_final=2.5,
            diagnostic_stride=500
        )
        solver = CoupledPhaseFieldSolver(cfg)
        res = solver.run(verbose=False)
        
        results[beta] = {
            'time': np.array(res['history']['time']),
            'dissolution': np.array(res['history']['dissolution_fraction']),
            'release': np.array(res['history']['release_fraction']),
            'radius': np.array(res['history']['effective_radius'])
        }
        
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(beta_values)))
    
    for i, beta in enumerate(beta_values):
        t = results[beta]['time']
        diss = results[beta]['dissolution'] * 100.0
        rad = results[beta]['radius']
        ax1.plot(t, diss, label=fr'$\beta_{{shear}} = {beta:.2f}$', color=colors[i], linewidth=2.2)
        ax2.plot(t, rad, label=fr'$\beta_{{shear}} = {beta:.2f}$', color=colors[i], linewidth=2.2)
        
    ax1.set_title('A. Matrix Dissolution vs Fluid Shear Coupling', fontsize=11, fontweight='bold')
    ax1.set_xlabel('Dimensionless Time (t*)', fontsize=10)
    ax1.set_ylabel('Lipid Matrix Degraded (%)', fontsize=10)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(frameon=True)
    
    ax2.set_title(r'B. Effective Nanocarrier Radius $R_{eff}(t)$', fontsize=11, fontweight='bold')
    ax2.set_xlabel('Dimensionless Time (t*)', fontsize=10)
    ax2.set_ylabel(r'Core Radius $R_{eff} / R_0$', fontsize=10)
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend(frameon=True)
    
    plt.suptitle('Hydrodynamic Shear Stress Coupling and Nanocarrier Degradation', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{save_dir}/figures/parametric_shear_coupling.png", dpi=300)
    plt.close()
    
    print(f"Completed Shear sweep. Figure saved to {save_dir}/figures/parametric_shear_coupling.png")
    return results


def run_cahn_number_sweep(cn_values=[0.18, 0.25, 0.35], save_dir="results"):
    """Evaluates the effect of diffuse interface width (Cahn number Cn) on release."""
    print("\n========================================================")
    print("  PARAMETRIC SWEEP 3: CAHN NUMBER (Cn = epsilon / R0) VARIATION")
    print("========================================================")
    
    os.makedirs(f"{save_dir}/figures", exist_ok=True)
    results = {}
    
    for cn in cn_values:
        print(f"--> Simulating Cn = {cn:4.2f} ...")
        cfg = SimulationConfig(
            Lx=20.0, Ly=10.0, Nx=100, Ny=50,
            R0=2.2, xc=5.0, yc=5.0,
            u_max=1.5,
            Cn=cn,
            k_diss=0.08,
            Pe_c=30.0,
            dt=2e-4,
            t_final=2.0,
            diagnostic_stride=500
        )
        solver = CoupledPhaseFieldSolver(cfg)
        res = solver.run(verbose=False)
        
        results[cn] = {
            'time': np.array(res['history']['time']),
            'release': np.array(res['history']['release_fraction']),
            'radius': np.array(res['history']['effective_radius'])
        }
        
    plt.figure(figsize=(9, 5))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    for i, cn in enumerate(cn_values):
        t = results[cn]['time']
        rel = results[cn]['release'] * 100.0
        plt.plot(t, rel, label=f'Cn = {cn:.2f} (Interface width ~ {np.sqrt(2)*cn:.2f} R0)', color=colors[i], linewidth=2.2)
        
    plt.title('Effect of Interfacial Diffuseness (Cahn Number Cn) on Release Kinetics', fontsize=12, fontweight='bold')
    plt.xlabel('Dimensionless Time (t*)', fontsize=11)
    plt.ylabel('Fractional Drug Released (%)', fontsize=11)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(frameon=True, fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{save_dir}/figures/parametric_cahn_number.png", dpi=300)
    plt.close()
    
    print(f"Completed Cahn number sweep. Figure saved to {save_dir}/figures/parametric_cahn_number.png")
    return results


if __name__ == "__main__":
    t0 = time.time()
    res_pe = run_peclet_sweep()
    res_shear = run_shear_sweep()
    res_cn = run_cahn_number_sweep()
    print(f"\nAll parametric sweeps successfully executed in {time.time()-t0:.2f} seconds!")

