"""
fit_M1_dose_response.py

Fits Hill-function copper dose-response parameters (K, n) from real M1 plate-reader
data and compares against the literature-proxy baseline used in Module 1
(K_COPPER = 5.0 uM, N_HILL = 1.0, from Fu et al. 2024 / CueR-PcopA context).

Usage:
    python fit_M1_dose_response.py M1_dose_response_data.csv

Input CSV columns (see M1_dose_response_template.csv):
    construct, replicate_id, cu_uM, time_h, OD600, GFP_raw, GFP_OD_normalized,
    operator, date, plate_position

Output:
    - m1_fit_results.csv     : fitted K, n, beta_max, beta_leak with 95% CI
    - m1_fit_comparison.png  : real data vs literature-proxy curve overlay
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Literature-proxy baseline (Module 1 grounded parameters)
LIT_K = 5.0        # uM, CueR/P(copA) half-max
LIT_N = 1.0         # Hill coefficient
LIT_BETA_MAX = 1000.0
LIT_BETA_LEAK = 100.0   # 10% of beta_max, Fu et al. 2024


def hill(cu, beta_max, beta_leak, K, n):
    return beta_leak + (beta_max - beta_leak) * (cu**n) / (K**n + cu**n)


def fit_dose_response(df, time_h=None):
    """Fit Hill parameters to GFP_OD_normalized vs cu_uM at a chosen readout time.
    If time_h is None, uses the most common time_h in the data."""
    if time_h is None:
        time_h = df["time_h"].mode().iloc[0]
    sub = df[df["time_h"] == time_h].copy()
    if len(sub) < 4:
        raise ValueError(f"Need >=4 points at time_h={time_h} to fit 4 parameters; got {len(sub)}")
    x = sub["cu_uM"].values.astype(float)
    y = sub["GFP_OD_normalized"].values.astype(float)

    p0 = [y.max(), y.min(), LIT_K, LIT_N]
    bounds = ([0, 0, 1e-3, 0.3], [np.inf, np.inf, 1000, 4])
    popt, pcov = curve_fit(hill, x, y, p0=p0, bounds=bounds, maxfev=10000)
    perr = np.sqrt(np.diag(pcov))
    beta_max, beta_leak, K, n = popt
    return {
        "time_h": time_h, "n_points": len(sub),
        "beta_max": beta_max, "beta_max_se": perr[0],
        "beta_leak": beta_leak, "beta_leak_se": perr[1],
        "K_uM": K, "K_se": perr[2],
        "n_hill": n, "n_se": perr[3],
    }, x, y


def make_comparison_plot(fit, x, y, outpath="m1_fit_comparison.png"):
    cu_grid = np.logspace(-2, 3, 200)
    y_fit = hill(cu_grid, fit["beta_max"], fit["beta_leak"], fit["K_uM"], fit["n_hill"])
    y_lit = hill(cu_grid, LIT_BETA_MAX, LIT_BETA_LEAK, LIT_K, LIT_N)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(cu_grid, y_lit, "--", color="#888888", label=f"literature proxy (K={LIT_K}, n={LIT_N})")
    ax.plot(cu_grid, y_fit, "-", color="#1f6fb4",
            label=f"measured fit (K={fit['K_uM']:.2f}+/-{fit['K_se']:.2f}, n={fit['n_hill']:.2f}+/-{fit['n_se']:.2f})")
    ax.scatter(x, y, color="#1f6fb4", zorder=5, label="measured data")
    ax.set_xscale("log")
    ax.set_xlabel("copper (uM, log)")
    ax.set_ylabel("GFP/OD600 (normalized)")
    ax.set_title(f"M1 dose-response: measured vs literature-proxy fit (t={fit['time_h']}h)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    return fig


def main(csv_path):
    df = pd.read_csv(csv_path, comment="#")
    fit, x, y = fit_dose_response(df)
    pd.DataFrame([fit]).to_csv("m1_fit_results.csv", index=False)
    make_comparison_plot(fit, x, y)
    print(f"Fitted K={fit['K_uM']:.3f} uM (lit proxy {LIT_K}), n={fit['n_hill']:.3f} (lit proxy {LIT_N})")
    print("Saved m1_fit_results.csv and m1_fit_comparison.png")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "M1_dose_response_template.csv")
