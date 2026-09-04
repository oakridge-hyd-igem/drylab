"""
m1_uncertainty_band.py

Monte Carlo prediction band for the M1 copper dose-response curve, built from the
grounded parameter plausible ranges already defined in Module 7's sensitivity analysis
(PARAMS table: beta_leak0, K_copper, n_hill). This turns "measured data looks like it
matches the model" into a checkable claim: does the measured curve fall inside the
model's own 90% prediction band, given realistic parameter uncertainty?

This does NOT replace fit_M1_dose_response.py (which fits real data to get point
estimates of K and n). This script instead answers: "before we had real data, what
range of dose-response curves did the literature-grounded uncertainty already predict,
and did the measurement fall inside that range or outside it?"

Caveat: the band reflects PARAMETER uncertainty only (beta_leak0, K, n drawn from their
grounded plausible ranges), not measurement noise. Near saturation (high copper), nearly
all parameter draws converge to beta_max, so the band narrows to near-zero width there
even though real measurement noise does not shrink. A measured point falling outside the
band at high copper is not automatically a model failure; check whether it is within the
plate-reader's own noise floor first.

Usage:
    python m1_uncertainty_band.py                          # band only, no real data
    python m1_uncertainty_band.py M1_dose_response_data.csv # band + measured overlay + coverage check

Output:
    - m1_uncertainty_band.png   : median curve + 90% and 50% prediction bands
    - m1_band_coverage.csv      : (if real data given) fraction of measured points inside band
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RNG_SEED = 0
N_DRAWS = 5000

# Grounded plausible ranges, from module7_sensitivity.py PARAMS table (matches
# Modules 1-5 as saved). beta_max held fixed: it is the strong-RBS induced ceiling,
# a design choice, not a measured/uncertain biological parameter.
BETA_MAX = 1000.0
RANGES = {
    "beta_leak0": (10.0, 500.0),   # T1: Fu 2024 copper-promoter fold-change grounding
    "K_copper":   (1.0, 15.0),     # sensor half-max, inferred from operating windows
    "n_hill":     (1.0, 2.0),      # Stoyanov abstract + MerR-family -> ~1, swept to 2
}
BASELINE = {"beta_leak0": 100.0, "K_copper": 5.0, "n_hill": 1.0}


def hill(cu, beta_max, beta_leak, K, n):
    return beta_leak + (beta_max - beta_leak) * (cu**n) / (K**n + cu**n)


def monte_carlo_band(cu_grid, n_draws=N_DRAWS, seed=RNG_SEED):
    rng = np.random.default_rng(seed)
    draws = np.empty((n_draws, len(cu_grid)))
    bl_lo, bl_hi = RANGES["beta_leak0"]
    K_lo, K_hi = RANGES["K_copper"]
    n_lo, n_hi = RANGES["n_hill"]
    beta_leak_samples = rng.uniform(bl_lo, bl_hi, n_draws)
    K_samples = rng.uniform(K_lo, K_hi, n_draws)
    n_samples = rng.uniform(n_lo, n_hi, n_draws)
    for i in range(n_draws):
        draws[i] = hill(cu_grid, BETA_MAX, beta_leak_samples[i], K_samples[i], n_samples[i])
    return draws


def band_stats(draws):
    return {
        "p5": np.percentile(draws, 5, axis=0),
        "p25": np.percentile(draws, 25, axis=0),
        "p50": np.percentile(draws, 50, axis=0),
        "p75": np.percentile(draws, 75, axis=0),
        "p95": np.percentile(draws, 95, axis=0),
    }


def make_plot(cu_grid, stats, measured=None, outpath="m1_uncertainty_band.png"):
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.fill_between(cu_grid, stats["p5"], stats["p95"], alpha=0.2, color="#1f6fb4",
                     label="90% prediction band (parameter uncertainty)")
    ax.fill_between(cu_grid, stats["p25"], stats["p75"], alpha=0.35, color="#1f6fb4",
                     label="50% prediction band")
    ax.plot(cu_grid, stats["p50"], "-", color="#1f6fb4", label="median prediction")
    baseline_curve = hill(cu_grid, BETA_MAX, BASELINE["beta_leak0"], BASELINE["K_copper"], BASELINE["n_hill"])
    ax.plot(cu_grid, baseline_curve, "--", color="#888888", label="point-estimate baseline (K=5, n=1)")
    if measured is not None:
        ax.scatter(measured["cu_uM"], measured["GFP_OD_normalized"], color="#c0392b",
                   zorder=5, label="measured data", s=30)
    ax.set_xscale("log")
    ax.set_xlabel("copper (uM, log)")
    ax.set_ylabel("GFP/OD600 (normalized)")
    ax.set_title("M1 dose-response: literature-grounded prediction band vs measurement")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    return fig


def check_coverage(df, cu_grid, stats):
    """For each measured point, check whether it falls inside the 90% band at its cu_uM."""
    rows = []
    for _, row in df.iterrows():
        cu = row["cu_uM"]
        idx = np.argmin(np.abs(cu_grid - cu)) if cu > 0 else 0
        lo, hi = stats["p5"][idx], stats["p95"][idx]
        inside = lo <= row["GFP_OD_normalized"] <= hi
        rows.append({"cu_uM": cu, "measured": row["GFP_OD_normalized"],
                     "band_p5": lo, "band_p95": hi, "inside_90pct_band": inside})
    out = pd.DataFrame(rows)
    return out


def main(csv_path=None):
    cu_grid = np.logspace(-2, 3, 200)
    draws = monte_carlo_band(cu_grid)
    stats = band_stats(draws)

    measured = None
    if csv_path:
        df = pd.read_csv(csv_path, comment="#")
        measured = df
        coverage = check_coverage(df, cu_grid, stats)
        coverage.to_csv("m1_band_coverage.csv", index=False)
        frac_inside = coverage["inside_90pct_band"].mean()
        print(f"Coverage: {frac_inside:.1%} of measured points fall inside the 90% prediction band")
        print("Saved m1_band_coverage.csv")

    make_plot(cu_grid, stats, measured)
    print("Saved m1_uncertainty_band.png")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
