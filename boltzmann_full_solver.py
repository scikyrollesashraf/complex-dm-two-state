#!/usr/bin/env python3
# boltzmann_full_solver.py
# Scientific Boltzmann solver + two-state P(t) + parameter scan
# Author: prepared for Kyrolles Ashraf Dawood
# Requires: numpy, scipy, matplotlib

import os
import numpy as np
import math
from math import pi
from scipy import special
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import csv

# -----------------------------
# Physical constants (GeV-based units)
# -----------------------------
# NOTE: all internal energy units are GeV, length units GeV^{-1}, time units GeV^{-1}
hbar_GeV_s = 6.582119569e-25      # ℏ in GeV*s  (for converting to seconds when needed)
Mpl = 1.2209e19                   # Planck mass in GeV (not reduced)
GeV_to_kg = 1.78266192e-27        # not used directly here, but provided for reference

# Conversion: 1 cm^3 / s  -> GeV^{-2}
# Derived from: 1 cm = 5.067730652e13 GeV^{-1}; 1 s = 1.51926744e24 GeV^{-1}
cm_to_GeVinv = 5.067730652e13
s_to_GeVinv = 1.519267444e24
cm3s_to_GeV2 = (cm_to_GeVinv**3) / s_to_GeVinv
# numeric value ~ 8.57e16
#print("cm3s_to_GeV2 =", cm3s_to_GeV2)

# Numeric conversion factor for Omega h^2 from yield:
# Omega h^2 = (m [GeV]) * s0 [cm^-3] * Y_inf / (rho_crit/h^2 [GeV cm^-3])
# We use the standard numeric: 2.744e8 * (m/GeV) * Y
Omega_prefactor = 2.744e8

# -----------------------------
# Approximate g_star(T) and g_star_s(T) tables
# (Important: these are approximate Standard Model values; for publication use a full table)
# -----------------------------
# Tabulated points (T in GeV) and approximate g*
T_table = np.array([
    1e-10, 1e-4, 1e-3, 1e-2, 1e-1,
    0.2, 0.5, 1.0, 3.0, 10.0, 100.0, 1000.0
])
# approximate g_star (effective degrees of freedom for energy) — piecewise, approximate
gstar_table = np.array([
    3.36,   # very low T (photons + neutrinos)
    3.36,
    3.36,
    10.75,  # neutrino decoupling region
    17.25,  # QCD transition region (approx)
    61.75,
    75.0,
    86.25,
    95.0,
    96.25,
    106.75,
    106.75
])
# take gstar_s ~ gstar for simplicity in this code (approx)
gstar_s_table = gstar_table.copy()

def gstar(T):
    """Return approximate g_* (energy) for temperature T in GeV (linear interpolation in log-T)."""
    # avoid zero or negative T
    if T <= T_table[0]:
        return float(gstar_table[0])
    if T >= T_table[-1]:
        return float(gstar_table[-1])
    # linear interpolation in log T
    logT = np.log10(T)
    log_table = np.log10(T_table)
    return float(10**np.interp(logT, log_table, np.log10(gstar_table)) if False else np.interp(logT, log_table, gstar_table))

def gstar_s(T):
    """Return approximate g_*s (entropy) for temperature T in GeV."""
    if T <= T_table[0]:
        return float(gstar_s_table[0])
    if T >= T_table[-1]:
        return float(gstar_s_table[-1])
    logT = np.log10(T)
    log_table = np.log10(T_table)
    return float(np.interp(logT, log_table, gstar_s_table))

# -----------------------------
# Thermodynamic functions
# -----------------------------
def s_entropy(T):
    """Entropy density s(T) in units GeV^3 (natural units)"""
    g_s = gstar_s(T)
    # s = (2*pi^2/45) * g_{*s} * T^3
    return (2.0 * pi**2 / 45.0) * g_s * T**3

def Hubble_rad(T):
    """Hubble rate in radiation era, H(T) in GeV (natural units: 1 GeV = 1/(1.52e-24 s))"""
    g = gstar(T)
    # H = 1.66 * sqrt(g_*) * T^2 / Mpl   (Mpl in GeV)
    return 1.66 * math.sqrt(g) * T**2 / Mpl

# -----------------------------
# Equilibrium number density and equilibrium yield
# -----------------------------
def n_eq_nonrel(m, T, g=1):
    """Number density n_eq for a nonrelativistic Maxwell-Boltzmann particle:
       n_eq = g * (m^2 T / (2*pi^2)) * K2(m/T)
       returns n_eq in units GeV^3 (natural units)."""
    x = m / T
    # use Bessel K2 from scipy.special
    K2 = special.kv(2, x)
    pref = g * (m**2 * T) / (2.0 * pi**2)
    return pref * K2

def Y_eq(m, T, g=1):
    """Equilibrium yield Y_eq = n_eq / s"""
    s = s_entropy(T)
    if s <= 0:
        return 0.0
    return n_eq_nonrel(m, T, g) / s

# -----------------------------
# Boltzmann RHS
# -----------------------------
def boltzmann_rhs(x, Y, params):
    """
    RHS for dY/dx.
    x = m/T
    Y: [Y] (array-like, but we use scalar)
    params: dict containing keys:
      - m (GeV), sigma_v (GeV^-2), Gamma_decay (GeV), S_fi (function or scalar), g (internal dof)
    returns dY/dx (scalar)
    """
    m = params['m']
    sigma_v = params['sigma_v']   # in GeV^-2
    Gamma = params.get('Gamma_decay', 0.0)  # in GeV
    S_fi = params.get('S_fi', 0.0)           # source term (GeV^4?) expected as yield-rate term; here scalar in same units as Gamma for simplicity
    g = params.get('g', 1)
    T = m / x
    # thermodynamic quantities
    s = s_entropy(T)
    H = Hubble_rad(T)
    Yeq = Y_eq(m, T, g)
    # Avoid division by zero
    if H <= 0 or x <= 0:
        return 0.0
    # Boltzmann form:
    dYdx = - (s * sigma_v) / (x * H) * (Y[0]**2 - Yeq**2) - (Gamma / (x * H)) * Y[0] + (S_fi / (x * H))
    return dYdx

# -----------------------------
# Solver wrapper
# -----------------------------
def solve_Y(m, sigma_v_cm3s, Gamma_s_inv=0.0, S_fi=0.0, g=1, x_init=1.0, x_final=1e4, rtol=1e-8, atol=1e-12):
    """
    Solve Boltzmann for given parameters and return x array and Y array.
    sigma_v_cm3s: input in cm^3/s (commonly used); converted to GeV^-2 internally.
    Gamma_s_inv: decay width in s^-1 (optional) -- converted to GeV
    S_fi: freeze-in source term (interpreted here in GeV units scaled consistently) -- can be 0
    """
    # Unit conversions:
    sigma_v = sigma_v_cm3s * cm3s_to_GeV2   # convert to GeV^-2
    # convert decay width in s^-1 to GeV: 1 s^-1 = 1 / (1.519267444e24) GeV
    Gamma_GeV = Gamma_s_inv / s_to_GeVinv if Gamma_s_inv != 0.0 else 0.0
    params = {'m': m, 'sigma_v': sigma_v, 'Gamma_decay': Gamma_GeV, 'S_fi': S_fi, 'g': g}
    # initial condition Y(x_init) ~ Y_eq(x_init) for equilibrium start
    T_init = m / x_init
    Y0 = Y_eq(m, T_init, g)
    # solve
    def rhs(x, Y): 
        return boltzmann_rhs(x, Y, params)
    sol = solve_ivp(rhs, (x_init, x_final), [Y0], method='BDF', dense_output=True, rtol=rtol, atol=atol, max_step=1e2)
    x_arr = sol.t
    Y_arr = sol.y[0]
    return x_arr, Y_arr, params

# -----------------------------
# P(t) two-state plotting function
# -----------------------------
def P_r_to_i_time(t_s, V_eV, Delta_eV):
    """
    Compute P_{R->I}(t) for given V (eV) and Delta (eV).
    t_s: times array in seconds
    Note: convert eV -> GeV inside (1 eV = 1e-9 GeV)
    """
    # convert eV to GeV
    V = V_eV * 1e-9
    Delta = Delta_eV * 1e-9
    Omega = np.sqrt(Delta**2 + V**2)
    # convert t in s to natural units (GeV^-1): multiply by (1 s) = 1.519267444e24 GeV^-1
    t_GeVinv = np.array(t_s) * s_to_GeVinv
    arg = (Omega * t_GeVinv)  # dimensionless (Omega [GeV] * t [GeV^-1])
    P = (V**2 / (Delta**2 + V**2)) * np.sin(arg)**2
    return P

# -----------------------------
# Utilities for plotting and saving results
# -----------------------------
def plot_Pt_example(outdir='outputs', V_eV=1e-4, Delta_eV=0.0):
    os.makedirs(outdir, exist_ok=True)
    t = np.logspace(-12, 2, 2000)  # from 1e-12 s to 100 s (example)
    P = P_r_to_i_time(t, V_eV, Delta_eV)
    plt.figure(figsize=(6,4))
    plt.semilogx(t, P)
    plt.xlabel('t (s)')
    plt.ylabel(r'$P_{R\to I}(t)$')
    plt.title(f'P(t): V={V_eV} eV, Delta={Delta_eV} eV')
    plt.grid(True)
    fname = os.path.join(outdir, 'P_t.pdf')
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
    print("Saved P(t) to", fname)

def plot_Y_vs_x(x, Y, outdir='outputs'):
    os.makedirs(outdir, exist_ok=True)
    plt.figure(figsize=(6,4))
    plt.loglog(x, Y, label='Y(x)')
    plt.xlabel('x = m/T')
    plt.ylabel('Y (yield)')
    plt.title('Boltzmann solution Y(x)')
    plt.grid(True)
    plt.legend()
    fname = os.path.join(outdir, 'Y_vs_x.pdf')
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
    print("Saved Y(x) to", fname)

# -----------------------------
# Example / main runner
# -----------------------------
def run_demo():
    outdir = 'outputs'
    os.makedirs(outdir, exist_ok=True)
    # Example parameters
    m_GeV = 10.0                 # DM mass in GeV
    sigma_v_cm3s = 3e-26         # canonical thermal cross-section (cm^3/s)
    Gamma_s_inv = 0.0            # no decay
    # solve
    x, Y, params = solve_Y(m_GeV, sigma_v_cm3s, Gamma_s_inv, S_fi=0.0, g=1, x_init=1.0, x_final=1e4)
    Y_inf = float(Y[-1])
    Omega_h2 = Omega_prefactor * (m_GeV) * Y_inf
    print(f"Demo: m={m_GeV} GeV, sigma_v={sigma_v_cm3s} cm^3/s -> Y_inf={Y_inf:.3e}, Omega h^2 = {Omega_h2:.3e}")
    # save CSV
    csvfile = os.path.join(outdir, 'Y_vs_x.csv')
    with open(csvfile, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['x','Y'])
        for xi, yi in zip(x, Y):
            writer.writerow([xi, yi])
    print("Saved CSV:", csvfile)
    # plots
    plot_Y_vs_x(x, Y, outdir)
    plot_Pt_example(outdir, V_eV=1e-6, Delta_eV=0.0)
    # write summary
    with open(os.path.join(outdir, 'summary.txt'), 'w') as f:
        f.write(f"m={m_GeV} GeV\n")
        f.write(f"sigma_v (cm^3/s) = {sigma_v_cm3s}\n")
        f.write(f"Y_inf = {Y_inf:.6e}\n")
        f.write(f"Omega h^2 = {Omega_h2:.6e}\n")
    print("Demo complete. Outputs in", outdir)

if __name__ == '__main__':
    run_demo()
