"""
Pharmacokinetics, Dissolution Kinetics & Mass Transfer Analysis Suite
Performs:
1. Classical Pharmaceutical Model Fitting:
   - Higuchi Matrix Diffusion Model (F = K_H * t^0.5)
   - Korsmeyer-Peppas Anomalous Transport (ln F = ln K_KP + n * ln t)
   - Hixson-Crowell Shrinking Core Model ((1 - F)^(1/3) = 1 - K_HC * t)
   - First-Order Dissolution Model (F = 1 - exp(-k_1 * t))
2. Convective Mass Transfer & Sherwood Number Scaling:
   - Evaluates Sh vs Pe across low-to-moderate Reynolds creeping shear flow
   - Compares with theoretical Levich boundary-layer scaling (Sh ~ Pe^(1/3))
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import csv
from scipy.optimize import curve_fit

from src.coupled_solver import SimulationConfig, CoupledPhaseFieldSolver


def fit_pharmaceutical_models(t: np.ndarray, F: np.ndarray, save_dir: str = "results"):
    """
    Fits empirical and mechanistic release models to fractional drug release data F(t).
    Only fits up to 60-80% release where standard models are valid.
    """
    os.makedirs(f"{save_dir}/figures", exist_ok=True)
    os.makedirs(f"{save_dir}/data", exist_ok=True)
    
    # Filter valid range for fitting (0.02 < F < 0.85)
    valid_mask = (F >= 0.02) & (F <= 0.85) & (t > 0)
    t_fit = t[valid_mask]
    F_fit = F[valid_mask]
    
    if len(t_fit) < 4:
        t_fit = t[1:]
        F_fit = F[1:]

    # 1. Higuchi Model: F = K_H * sqrt(t)
    def higuchi_func(t_val, k_h):
        return k_h * np.sqrt(t_val)

    popt_h, _ = curve_fit(higuchi_func, t_fit, F_fit, p0=[0.4])
    k_h = popt_h[0]
    F_pred_h = higuchi_func(t_fit, k_h)
    ss_res_h = np.sum((F_fit - F_pred_h)**2)
    ss_tot = np.sum((F_fit - np.mean(F_fit))**2)
    r2_h = 1.0 - (ss_res_h / (ss_tot + 1e-12))
    
    # 2. Korsmeyer-Peppas Model: ln(F) = ln(K_kp) + n * ln(t) -> F = K_kp * t^n
    log_t = np.log(t_fit)
    log_F = np.log(F_fit)
    coeffs_kp = np.polyfit(log_t, log_F, 1)
    n_kp = coeffs_kp[0]
    k_kp = np.exp(coeffs_kp[1])
    F_pred_kp = k_kp * (t_fit ** n_kp)
    ss_res_kp = np.sum((F_fit - F_pred_kp)**2)
    r2_kp = 1.0 - (ss_res_kp / (ss_tot + 1e-12))
    
    # 3. Hixson-Crowell Model: (1 - F)^(1/3) = 1 - K_hc * t
    y_hc = 1.0 - (1.0 - F_fit)**(1.0/3.0)
    coeffs_hc = np.polyfit(t_fit, y_hc, 1)
    k_hc = coeffs_hc[0]
    F_pred_hc = 1.0 - (1.0 - k_hc * t_fit)**3.0
    ss_res_hc = np.sum((F_fit - F_pred_hc)**2)
    r2_hc = 1.0 - (ss_res_hc / (ss_tot + 1e-12))
    
    # 4. First-Order Model: F = 1 - exp(-k1 * t)
    def first_order_func(t_val, k1):
        return 1.0 - np.exp(-k1 * t_val)
    
    popt_fo, _ = curve_fit(first_order_func, t_fit, F_fit, p0=[0.5])
    k_fo = popt_fo[0]
    F_pred_fo = first_order_func(t_fit, k_fo)
    ss_res_fo = np.sum((F_fit - F_pred_fo)**2)
    r2_fo = 1.0 - (ss_res_fo / (ss_tot + 1e-12))
    
    # Export fitting results to CSV
    csv_path = f"{save_dir}/data/kinetics_model_summary.csv"
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Model", "Equation", "Parameters", "R_squared", "Transport Mechanism"])
        writer.writerow(["Higuchi", "F = K_H * t^(0.5)", f"K_H = {k_h:.4f}", f"{r2_h:.4f}", "Pure Fickian diffusion from matrix"])
        writer.writerow(["Korsmeyer-Peppas", "F = K_KP * t^n", f"K_KP = {k_kp:.4f}, n = {n_kp:.4f}", f"{r2_kp:.4f}", 
                         f"Anomalous transport ({'Case II erosion' if n_kp > 0.85 else 'Coupled diffusion-relaxation'})"])
        writer.writerow(["Hixson-Crowell", "(1-F)^(1/3) = 1 - K_HC * t", f"K_HC = {k_hc:.4f}", f"{r2_hc:.4f}", "Dissolution from shrinking sphere diameter"])
        writer.writerow(["First-Order", "F = 1 - exp(-k1 * t)", f"k1 = {k_fo:.4f}", f"{r2_fo:.4f}", "Concentration-gradient driven release"])
    print(f"Exported model fitting summary to {csv_path}")
    
    # Plotting comparison
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 9))
    
    # Ax1: Overall comparison
    t_dense = np.linspace(t_fit[0], t_fit[-1], 200)
    ax1.plot(t, F * 100.0, 'ko', markersize=4, label='Phase-Field Simulation Data', alpha=0.6)
    ax1.plot(t_dense, higuchi_func(t_dense, k_h) * 100.0, 'r--', linewidth=2.0, label=f'Higuchi (R² = {r2_h:.3f})')
    ax1.plot(t_dense, (k_kp * t_dense**n_kp) * 100.0, 'b-', linewidth=2.0, label=f'Korsmeyer-Peppas (R² = {r2_kp:.3f}, n = {n_kp:.2f})')
    ax1.plot(t_dense, (1.0 - (1.0 - k_hc * t_dense)**3) * 100.0, 'g-.', linewidth=2.0, label=f'Hixson-Crowell (R² = {r2_hc:.3f})')
    ax1.set_title('A. Pharmaceutical Kinetics Model Overlay', fontsize=11, fontweight='bold')
    ax1.set_xlabel('Dimensionless Time (t*)')
    ax1.set_ylabel('Cumulative Drug Released (%)')
    ax1.legend(frameon=True, fontsize=9)
    ax1.grid(True, linestyle='--', alpha=0.5)
    
    # Ax2: Higuchi linearization (F vs sqrt(t))
    ax2.plot(np.sqrt(t_fit), F_fit, 'rs', markersize=5, label='Simulation')
    ax2.plot(np.sqrt(t_dense), higuchi_func(t_dense, k_h), 'k-', label=f'Linear Fit: slope={k_h:.3f}')
    ax2.set_title(r'B. Higuchi Matrix Diffusion Model ($F$ vs $\sqrt{t^*}$)', fontsize=11, fontweight='bold')
    ax2.set_xlabel(r'Square Root of Time $\sqrt{t^*}$')
    ax2.set_ylabel('Fraction Released F')
    ax2.legend(frameon=True, fontsize=9)
    ax2.grid(True, linestyle='--', alpha=0.5)
    
    # Ax3: Korsmeyer-Peppas log-log linearization
    ax3.plot(log_t, log_F, 'bo', markersize=5, label='Simulation')
    ax3.plot(log_t, coeffs_kp[0]*log_t + coeffs_kp[1], 'k-', label=fr'Slope $n = {n_kp:.3f}$')
    ax3.set_title(r'C. Korsmeyer-Peppas Log-Log Plot ($\ln F$ vs $\ln t^*$)', fontsize=11, fontweight='bold')
    ax3.set_xlabel(r'$\ln(t^*)$')
    ax3.set_ylabel(r'$\ln(F)$')
    ax3.legend(frameon=True, fontsize=9)
    ax3.grid(True, linestyle='--', alpha=0.5)
    
    # Ax4: Hixson-Crowell cube root linearization
    ax4.plot(t_fit, (1.0 - F_fit)**(1.0/3.0), 'g^', markersize=5, label='Simulation')
    ax4.plot(t_dense, 1.0 - k_hc * t_dense, 'k-', label=fr'Slope $K_{{HC}} = {k_hc:.3f}$')
    ax4.set_title(r'D. Hixson-Crowell Model ($(1-F)^{1/3}$ vs $t^*$)', fontsize=11, fontweight='bold')
    ax4.set_xlabel('Dimensionless Time (t*)')
    ax4.set_ylabel(r'$(1 - F)^{1/3}$')
    ax4.legend(frameon=True, fontsize=9)
    ax4.grid(True, linestyle='--', alpha=0.5)
    
    plt.suptitle('Pharmacokinetic & Mathematical Dissolution Model Evaluation', fontsize=13, fontweight='bold')
    plt.tight_layout()
    fig_path = f"{save_dir}/figures/pharmacokinetics_models.png"
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"Saved kinetic model fitting plot to {fig_path}")
    
    return {
        'higuchi': {'k': k_h, 'r2': r2_h},
        'korsmeyer_peppas': {'k': k_kp, 'n': n_kp, 'r2': r2_kp},
        'hixson_crowell': {'k': k_hc, 'r2': r2_hc},
        'first_order': {'k': k_fo, 'r2': r2_fo}
    }


def analyze_sherwood_scaling(save_dir: str = "results"):
    """
    Evaluates Sherwood number Sh as a function of Péclet number Pe.
    Tests agreement with Levich boundary-layer theory (Sh ~ Pe^(1/3)).
    """
    print("\n========================================================")
    print("  MASS TRANSFER & SHERWOOD NUMBER (Sh vs Pe) SCALING")
    print("========================================================")
    
    pe_list = [5.0, 15.0, 35.0, 80.0, 200.0, 500.0]
    sh_list = []
    
    for pe in pe_list:
        print(f"--> Computing asymptotic Sh for Pe = {pe:5.1f} ...")
        cfg = SimulationConfig(
            Lx=20.0, Ly=10.0, Nx=100, Ny=50,
            R0=2.2, xc=5.0, yc=5.0,
            u_max=1.5,
            k_diss=0.08,
            Pe_c=pe,
            dt=2e-4,
            t_final=1.5,
            diagnostic_stride=500
        )
        solver = CoupledPhaseFieldSolver(cfg)
        res = solver.run(verbose=False)
        
        # Quasi-steady Sherwood number in middle release window (t between 0.5 and 1.2)
        sh_arr = np.array(res['history']['sherwood_number'])
        t_arr = np.array(res['history']['time'])
        mid_mask = (t_arr >= 0.4) & (t_arr <= 1.2)
        mean_sh = np.mean(sh_arr[mid_mask]) if np.any(mid_mask) else sh_arr[-1]
        sh_list.append(mean_sh)
        
    pe_arr = np.array(pe_list)
    sh_arr = np.array(sh_list)
    
    # Power-law fit: ln(Sh) = alpha * ln(Pe) + ln(C0)
    log_pe = np.log(pe_arr)
    log_sh = np.log(sh_arr)
    poly = np.polyfit(log_pe, log_sh, 1)
    alpha = poly[0]
    c0 = np.exp(poly[1])
    
    print(f"Fitted Power-Law Exponent: Sh = {c0:.3f} * Pe^({alpha:.3f})")
    print(f"Theoretical Levich Creeping Shear Flow Exponent: alpha = 0.333 (1/3)")
    
    # Save figure
    os.makedirs(f"{save_dir}/figures", exist_ok=True)
    plt.figure(figsize=(9, 6))
    plt.loglog(pe_arr, sh_arr, 'ro', markersize=8, label='Phase-Field Simulation Data')
    
    pe_dense = np.logspace(np.log10(pe_arr[0]), np.log10(pe_arr[-1]), 100)
    plt.loglog(pe_dense, c0 * (pe_dense ** alpha), 'k-', linewidth=2.0, 
               label=fr'Fitted Scaling: $Sh \approx {c0:.2f} \, Pe^{{{alpha:.3f}}}$')
    
    # Theoretical 1/3 slope comparison line
    c_lev = sh_arr[2] / (pe_arr[2]**(1.0/3.0))
    plt.loglog(pe_dense, c_lev * (pe_dense ** (1.0/3.0)), 'b--', linewidth=1.8, 
               label=r'Levich Boundary-Layer Asymptote: $Sh \propto Pe^{1/3}$')
    
    plt.title('Sherwood Number vs Péclet Number: Convective Boundary Layer Scaling', fontsize=12, fontweight='bold')
    plt.xlabel('Solute Péclet Number (Pe = U_max R0 / D_plasma)', fontsize=11)
    plt.ylabel('Sherwood Number (Sh = k_m * 2R0 / D_plasma)', fontsize=11)
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.legend(frameon=True, fontsize=10)
    plt.tight_layout()
    
    fig_path = f"{save_dir}/figures/sherwood_vs_peclet.png"
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"Saved Sherwood vs Péclet scaling plot to {fig_path}")
    
    return {'pe': pe_list, 'sh': sh_list, 'alpha': alpha, 'c0': c0}


if __name__ == "__main__":
    # Generate baseline data for model fitting
    print("Generating baseline run for pharmacokinetic model fitting...")
    cfg = SimulationConfig(t_final=3.0, dt=2e-4, Nx=100, Ny=50, k_diss=0.08, Pe_c=30.0, diagnostic_stride=500)
    solver = CoupledPhaseFieldSolver(cfg)
    res = solver.run(verbose=False)
    
    t = np.array(res['history']['time'])
    F = np.array(res['history']['release_fraction'])
    
    fit_pharmaceutical_models(t, F)
    analyze_sherwood_scaling()

