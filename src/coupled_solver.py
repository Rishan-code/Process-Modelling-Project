"""
Coupled 2D Phase-Field and Convective Mass Transfer Solver
Simulates liposomal nanocarrier dissolution and shear-mediated drug release
under laminar microvascular blood flow (Hagen-Poiseuille regime).

Physics resolved:
1. Convective Cahn-Hilliard equation for moving lipid-plasma interface (phi)
2. Interfacial dissolution reaction coupled to fluid wall shear stress (tau_w)
3. Convective diffusion equation for encapsulated drug (C)
4. Phase-dependent mobility and effective diffusivity
5. Rigorous mass conservation: Carrier core + Plasma plume + Cumulative outflow
"""

import numpy as np
import matplotlib.pyplot as plt
import time
from dataclasses import dataclass
from typing import Dict, Any, Tuple, Optional


@dataclass
class SimulationConfig:
    # Spatial domain
    Lx: float = 24.0          # Channel length (normalized by R0)
    Ly: float = 12.0          # Channel height
    Nx: int = 120             # Grid points in x
    Ny: int = 60              # Grid points in y
    
    # Nanocarrier geometry
    R0: float = 2.5           # Initial liposome radius
    xc: float = 6.0           # Liposome initial center x
    yc: float = 6.0           # Liposome initial center y (channel midline)
    
    # Dimensionless phase-field parameters
    Cn: float = 0.25          # Cahn number (epsilon / R0)
    W: float = 1.0            # Double-well potential barrier height
    M0: float = 0.5           # Baseline mobility constant
    M_reg: float = 1e-4       # Regularization mobility
    
    # Flow & Dissolution
    u_max: float = 1.5        # Centerline flow velocity
    k_diss: float = 0.08      # Interfacial dissolution kinetic rate constant
    beta_shear: float = 0.5   # Shear-coupling coefficient for dissolution acceleration
    advect_carrier: bool = False  # False = wall-adhered / targeted nanocarrier; True = freely suspended
    
    # Drug transport
    Pe_c: float = 50.0        # Péclet number for drug transport (U_max * R0 / D_plasma)
    D_ratio: float = 1e-2     # Diffusivity contrast: D_core / D_plasma
    
    # Temporal integration
    dt: float = 2e-4          # Time step size
    t_final: float = 5.0      # Total non-dimensional simulation time
    diagnostic_stride: int = 500  # Steps between metric records


class CoupledPhaseFieldSolver:
    def __init__(self, config: SimulationConfig):
        self.cfg = config
        
        # Grid setup
        self.dx = self.cfg.Lx / (self.cfg.Nx - 1)
        self.dy = self.cfg.Ly / (self.cfg.Ny - 1)
        self.x = np.linspace(0, self.cfg.Lx, self.cfg.Nx)
        self.y = np.linspace(0, self.cfg.Ly, self.cfg.Ny)
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='ij')
        
        # Initialize velocity field (Hagen-Poiseuille flow)
        # u_x(y) = 4 * u_max * (y / Ly) * (1 - y / Ly)
        self.vx = 4.0 * self.cfg.u_max * (self.Y / self.cfg.Ly) * (1.0 - self.Y / self.cfg.Ly)
        self.vy = np.zeros_like(self.vx)
        
        # Shear rate magnitude |du/dy|
        self.shear_rate = (4.0 * self.cfg.u_max / self.cfg.Ly) * np.abs(1.0 - 2.0 * self.Y / self.cfg.Ly)
        max_shear = np.max(self.shear_rate) if np.max(self.shear_rate) > 0 else 1.0
        self.shear_norm = self.shear_rate / max_shear
        
        # Effective interface parameter
        self.epsilon = self.cfg.Cn * self.cfg.R0
        
        # Initialize fields
        self.reset()

    def reset(self):
        """Initializes phi, C, and tracking metrics."""
        cfg = self.cfg
        r = np.sqrt((self.X - cfg.xc)**2 + (self.Y - cfg.yc)**2)
        
        # Phase field: phi = 0 inside liposome, phi = 1 in blood plasma
        self.phi = 0.5 * (1.0 + np.tanh((r - cfg.R0) / (np.sqrt(2.0) * self.epsilon)))
        
        # Drug concentration field: C = 1 inside core, smoothly dropping to 0 outside
        self.C = np.clip(1.0 - self.phi, 0.0, 1.0)
        
        # Metrics history
        self.time = 0.0
        initial_lipid = float(np.sum(1.0 - self.phi) * self.dx * self.dy)
        initial_carrier_drug = float(np.sum((1.0 - self.phi) * self.C) * self.dx * self.dy)
        initial_total_drug = float(np.sum(self.C) * self.dx * self.dy)
        
        self.initial_lipid_mass = initial_lipid
        self.initial_carrier_drug_mass = initial_carrier_drug
        self.initial_total_drug = initial_total_drug
        
        self.history = {
            'time': [0.0],
            'lipid_mass': [initial_lipid],
            'carrier_drug': [initial_carrier_drug],
            'plasma_drug': [0.0],
            'cumulative_outflow': [0.0],
            'total_system_drug': [initial_total_drug],
            'effective_radius': [cfg.R0],
            'dissolution_fraction': [0.0],
            'release_fraction': [0.0],
            'sherwood_number': [0.0]
        }

    # Spatial differential operators (conservative 2nd-order FDM)
    def laplacian(self, field: np.ndarray) -> np.ndarray:
        """2nd-order central difference Laplacian with zero-flux Neumann BCs."""
        lap = np.zeros_like(field)
        dx2, dy2 = self.dx**2, self.dy**2
        
        lap[1:-1, 1:-1] = (
            (field[2:, 1:-1] - 2.0 * field[1:-1, 1:-1] + field[:-2, 1:-1]) / dx2 +
            (field[1:-1, 2:] - 2.0 * field[1:-1, 1:-1] + field[1:-1, :-2]) / dy2
        )
        
        # Neumann zero-flux BCs at boundaries
        lap[0, :] = lap[1, :]
        lap[-1, :] = lap[-2, :]
        lap[:, 0] = lap[:, 1]
        lap[:, -1] = lap[:, -2]
        return lap

    def grad(self, field: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """2nd-order central gradient with zero-gradient boundaries."""
        gx = np.zeros_like(field)
        gy = np.zeros_like(field)
        
        gx[1:-1, :] = (field[2:, :] - field[:-2, :]) / (2.0 * self.dx)
        gy[:, 1:-1] = (field[:, 2:] - field[:, :-2]) / (2.0 * self.dy)
        
        # Boundaries
        gx[0, :] = (field[1, :] - field[0, :]) / self.dx
        gx[-1, :] = (field[-1, :] - field[-2, :]) / self.dx
        gy[:, 0] = (field[:, 1] - field[:, 0]) / self.dy
        gy[:, -1] = (field[:, -1] - field[:, -2]) / self.dy
        return gx, gy

    def div_flux(self, jx: np.ndarray, jy: np.ndarray) -> np.ndarray:
        """Conservative divergence of vector flux (jx, jy)."""
        div = np.zeros_like(jx)
        div[1:-1, :] += (jx[2:, :] - jx[:-2, :]) / (2.0 * self.dx)
        div[:, 1:-1] += (jy[:, 2:] - jy[:, :-2]) / (2.0 * self.dy)
        
        # Boundary flux divergence
        div[0, :] = (jx[1, :] - jx[0, :]) / self.dx
        div[-1, :] = (jx[-1, :] - jx[-2, :]) / self.dx
        div[:, 0] = (jy[:, 1] - jy[:, 0]) / self.dy
        div[:, -1] = (jy[:, -1] - jy[:, -2]) / self.dy
        return div

    def advection_upwind(self, field: np.ndarray) -> np.ndarray:
        """First-order stable upwind advection term vx * d(field)/dx."""
        adv = np.zeros_like(field)
        adv[1:, :] = self.vx[1:, :] * (field[1:, :] - field[:-1, :]) / self.dx
        adv[0, :] = self.vx[0, :] * (field[1, :] - field[0, :]) / self.dx
        return adv

    def step(self) -> float:
        """Advances the coupled system by dt. Returns instantaneous outflow flux."""
        cfg = self.cfg
        dt = cfg.dt
        phi = self.phi
        C = self.C
        
        # 1. Cahn-Hilliard Phase-Field Evolution
        df_dphi = 2.0 * cfg.W * phi * (1.0 - phi) * (1.0 - 2.0 * phi)
        mu = df_dphi - (self.epsilon**2) * self.laplacian(phi)
        
        mob = cfg.M0 * (phi**2) * ((1.0 - phi)**2) + cfg.M_reg
        g_mu_x, g_mu_y = self.grad(mu)
        jx_phi = -mob * g_mu_x
        jy_phi = -mob * g_mu_y
        div_j_phi = self.div_flux(jx_phi, jy_phi)
        
        # Interfacial dissolution reaction: localized at diffuse interface phi*(1-phi)
        shear_factor = 1.0 + cfg.beta_shear * self.shear_norm
        r_diss = cfg.k_diss * phi * (1.0 - phi) * shear_factor
        
        adv_phi = self.advection_upwind(phi) if cfg.advect_carrier else 0.0
        dphi_dt = -div_j_phi - adv_phi + r_diss
        new_phi = np.clip(phi + dt * dphi_dt, 0.0, 1.0)
        
        # 2. Convective Diffusion of Encapsulated Drug
        h_phi = phi**2 * (3.0 - 2.0 * phi)
        D_eff = (cfg.D_ratio + (1.0 - cfg.D_ratio) * h_phi) / cfg.Pe_c
        
        g_c_x, g_c_y = self.grad(C)
        jx_c = -D_eff * g_c_x
        jy_c = -D_eff * g_c_y
        div_j_c = self.div_flux(jx_c, jy_c)
        adv_c = self.advection_upwind(C)
        
        dC_dt = -div_j_c - adv_c
        new_C = np.clip(C + dt * dC_dt, 0.0, 1.0)
        
        # 3. Boundary Conditions
        # Inflow (x = 0): clean incoming blood
        new_phi[0, :] = 1.0
        new_C[0, :] = 0.0
        
        # Outlet (x = Lx): convective Neumann outflow
        new_phi[-1, :] = new_phi[-2, :]
        new_C[-1, :] = new_C[-2, :]
        
        # Upper & lower walls: no-slip, zero flux
        new_phi[:, 0] = new_phi[:, 1]
        new_phi[:, -1] = new_phi[:, -2]
        new_C[:, 0] = new_C[:, 1]
        new_C[:, -1] = new_C[:, -2]
        
        # Instantaneous drug outflow rate at x = Lx
        outflow_rate = float(np.sum(self.vx[-1, :] * new_C[-1, :]) * self.dy)
        
        # Commit state
        self.phi = new_phi
        self.C = new_C
        self.time += dt
        
        return outflow_rate

    def run(self, verbose: bool = False) -> Dict[str, Any]:
        """Runs the simulation until t_final, logging metrics."""
        n_total_steps = int(np.ceil(self.cfg.t_final / self.cfg.dt))
        t_start = time.time()
        
        cum_outflow = self.history['cumulative_outflow'][-1]
        
        for step in range(1, n_total_steps + 1):
            outflow_rate = self.step()
            cum_outflow += outflow_rate * self.cfg.dt
            
            if step % self.cfg.diagnostic_stride == 0 or step == n_total_steps:
                current_lipid = float(np.sum(1.0 - self.phi) * self.dx * self.dy)
                current_carrier_drug = float(np.sum((1.0 - self.phi) * self.C) * self.dx * self.dy)
                current_plasma_drug = float(np.sum(self.phi * self.C) * self.dx * self.dy)
                total_system_drug = current_carrier_drug + current_plasma_drug + cum_outflow
                
                # Degradation metrics
                r_eff = float(np.sqrt(np.maximum(current_lipid / np.pi, 0.0)))
                diss_frac = float((self.initial_lipid_mass - current_lipid) / (self.initial_lipid_mass + 1e-12))
                rel_frac = float((self.initial_carrier_drug_mass - current_carrier_drug) / (self.initial_carrier_drug_mass + 1e-12))
                
                # Sherwood number: Sh = (dM_rel/dt * 2 * R0) / (2 * pi * R0 * DeltaC * D_plasma)
                dt_diag = self.time - self.history['time'][-1] + 1e-12
                rel_rate = (rel_frac - self.history['release_fraction'][-1]) / dt_diag
                sh_val = float(rel_rate * self.cfg.Pe_c * self.cfg.R0 / (2.0 * np.pi + 1e-6))
                
                self.history['time'].append(self.time)
                self.history['lipid_mass'].append(current_lipid)
                self.history['carrier_drug'].append(current_carrier_drug)
                self.history['plasma_drug'].append(current_plasma_drug)
                self.history['cumulative_outflow'].append(cum_outflow)
                self.history['total_system_drug'].append(total_system_drug)
                self.history['effective_radius'].append(r_eff)
                self.history['dissolution_fraction'].append(np.clip(diss_frac, 0.0, 1.0))
                self.history['release_fraction'].append(np.clip(rel_frac, 0.0, 1.0))
                self.history['sherwood_number'].append(sh_val)
                
                if verbose and (step % (self.cfg.diagnostic_stride * 5) == 0 or step == n_total_steps):
                    pct = 100.0 * step / n_total_steps
                    print(f"[{pct:5.1f}%] t={self.time:5.2f} | R_eff={r_eff:.3f} | Lipid Dissolved: {diss_frac*100:4.1f}% | Drug Released: {rel_frac*100:4.1f}%")
        
        elapsed = time.time() - t_start
        if verbose:
            print(f"Simulation completed in {elapsed:.2f} seconds ({n_total_steps} steps).")
            
        return {
            'history': self.history,
            'phi': self.phi,
            'C': self.C,
            'X': self.X,
            'Y': self.Y,
            'vx': self.vx,
            'vy': self.vy,
            'elapsed_sec': elapsed
        }

    def plot_summary(self, save_path: str = "results/figures/coupled_simulation_summary.png"):
        """Generates a comprehensive 4-panel publication-grade overview."""
        import os
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        
        # 1. Phase field & interface contour
        ax1 = axes[0, 0]
        c1 = ax1.pcolormesh(self.X, self.Y, self.phi, cmap='viridis', shading='gouraud', vmin=0, vmax=1)
        ax1.contour(self.X, self.Y, self.phi, levels=[0.5], colors='red', linewidths=2.0)
        plt.colorbar(c1, ax=ax1, label=r'Phase Field $\phi$ (0: Lipid, 1: Plasma)')
        ax1.set_title('A. Liposome Phase Field & Moving Interface Contour', fontsize=11, fontweight='bold')
        ax1.set_xlabel('Channel Length (x / R0)')
        ax1.set_ylabel('Channel Height (y / R0)')
        
        # 2. Convective Blood Flow & Streamlines
        ax2 = axes[0, 1]
        c2 = ax2.pcolormesh(self.X, self.Y, self.vx, cmap='plasma', shading='gouraud')
        # Overlay velocity streamlines
        stride = 2
        ax2.streamplot(self.x, self.y, self.vx.T, self.vy.T, color='white', density=0.8, linewidth=0.8, arrowsize=1.0)
        plt.colorbar(c2, ax=ax2, label=r'Velocity $u_x(y) / U_{max}$')
        ax2.set_title('B. Hagen-Poiseuille Convective Flow Field & Streamlines', fontsize=11, fontweight='bold')
        ax2.set_xlabel('Channel Length (x / R0)')
        ax2.set_ylabel('Channel Height (y / R0)')
        
        # 3. Drug Concentration Plume
        ax3 = axes[1, 0]
        c3 = ax3.pcolormesh(self.X, self.Y, self.C, cmap='magma', shading='gouraud', vmin=0, vmax=1)
        ax3.contour(self.X, self.Y, self.phi, levels=[0.5], colors='cyan', linewidths=1.5, linestyles='--')
        plt.colorbar(c3, ax=ax3, label=r'Drug Concentration $C / C_0$')
        ax3.set_title('C. Encapsulated Drug Plume & Downstream Boundary Layer', fontsize=11, fontweight='bold')
        ax3.set_xlabel('Channel Length (x / R0)')
        ax3.set_ylabel('Channel Height (y / R0)')
        
        # 4. Transient Release & Degradation Kinetics
        ax4 = axes[1, 1]
        t = np.array(self.history['time'])
        rel = np.array(self.history['release_fraction']) * 100.0
        diss = np.array(self.history['dissolution_fraction']) * 100.0
        rad = np.array(self.history['effective_radius']) / self.cfg.R0
        
        ax4.plot(t, rel, 'r-', linewidth=2.5, label='Drug Release Fraction (%)')
        ax4.plot(t, diss, 'b--', linewidth=2.0, label='Lipid Matrix Dissolution (%)')
        ax4.set_xlabel('Dimensionless Time (t*)', fontsize=11)
        ax4.set_ylabel('Fractional Conversion (%)', fontsize=11)
        ax4.set_ylim(-2, 102)
        ax4.grid(True, linestyle='--', alpha=0.6)
        
        ax4_twin = ax4.twinx()
        ax4_twin.plot(t, rad, 'g-.', linewidth=2.0, label=r'Normalized Radius $R_{eff} / R_0$')
        ax4_twin.set_ylabel(r'Core Radius $R_{eff} / R_0$', color='g', fontsize=11)
        ax4_twin.tick_params(axis='y', labelcolor='g')
        ax4_twin.set_ylim(0, 1.1)
        
        lines1, labels1 = ax4.get_legend_handles_labels()
        lines2, labels2 = ax4_twin.get_legend_handles_labels()
        ax4.legend(lines1 + lines2, labels1 + labels2, loc='center right', frameon=True)
        ax4.set_title('D. Transient Release & Degradation Kinetics', fontsize=11, fontweight='bold')
        
        plt.suptitle(f'Coupled Phase-Field Simulation of Liposomal Dissolution & Convective Drug Release\n(Pe = {self.cfg.Pe_c}, Cn = {self.cfg.Cn}, k_diss = {self.cfg.k_diss})', fontsize=13, fontweight='bold')
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"Summary figure saved successfully to {save_path}")
