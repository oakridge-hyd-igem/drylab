"""
fit_M2_pilot.py

Summarizes M2 pilot leak/switching data across glucose, no-inducer, and arabinose
conditions for both tag variants (untagged, LVA). Reports per the Rescue HQ Section 6
rule: ON-sequence detection frequency and OFF-only/mixed/ON-only categories, GFP output,
NOT a per-cell/per-molecule flip fraction.

Usage:
    python fit_M2_pilot.py M2_pilot_data.csv

Input CSV columns (see M2_pilot_template.csv):
    construct, tag_variant, induction, replicate_id, time_h, OD600, GFP_raw,
    GFP_OD_normalized, pcr_call, operator, date
    pcr_call in {OFF_only, mixed, ON_only}

Output:
    - m2_pilot_summary.csv    : per (tag_variant, induction, time_h) GFP stats + PCR category
                                fractions, plus the Dry Model D2 prediction and residual for
                                each condition's ON_only fraction (goodness-of-fit check)
    - m2_pilot_comparison.png : GFP/OD600 leak and switching, untagged vs LVA, by condition
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

VALID_CALLS = {"OFF_only", "mixed", "ON_only"}

# Dry Model D2 (re-mapped) predicted ON_only (flipped) fraction, from the grounded pBAD
# leak engine (D2_remap_model.py: k_flip=0.4/h, K_I=10, gamma_dil=0.3/h, gamma_LVA=1.04/h).
# Reproduced inline (not imported) so this template stays standalone.
_BETA_MAX_EFF = 1000.0 * 0.1  # strong-RBS ceiling * weak-RBS (B0033) factor
_K_FLIP, _K_I, _GAMMA_DIL, _GAMMA_LVA = 0.4, 10.0, 0.30, 1.04
_LEAK_FRAC = {"glucose": 1 / 1200, "no_inducer": 0.005, "arabinose": 1.0}


def _d2_predicted_flip(condition, tag_variant, t_h, k_flip=_K_FLIP, K_I=_K_I,
                        gamma_dil=_GAMMA_DIL, gamma_lva=_GAMMA_LVA):
    """Predicted ON_only (flipped) fraction at time t_h, from Dry Model D2."""
    gamma = gamma_dil + (gamma_lva if tag_variant == "LVA" else 0.0)
    beta = _LEAK_FRAC.get(condition, np.nan) * _BETA_MAX_EFF
    if np.isnan(beta):
        return np.nan
    dt = 1e-3
    tt = np.arange(0, t_h + dt, dt)
    I = (beta / gamma) * (1 - np.exp(-gamma * tt))
    lam = k_flip * I / (I + K_I)
    cum = np.concatenate([[0], np.cumsum((lam[1:] + lam[:-1]) / 2 * np.diff(tt))])
    return 1 - np.exp(-cum[-1])


# Plausible ranges for the flip-kinetics parameters, matching Module 7's PARAMS table
# (module7_sensitivity.py: k_flip, K_I, gamma_dil, gamma_LVA). The pBAD leak fractions
# themselves (glucose/no_inducer/arabinose) are held at their literature point estimates
# (Guzman 1995 / assumed) since no documented plausible range exists for them yet -
# only the flip-kinetics uncertainty is propagated here.
_KINETIC_RANGES = {
    "k_flip": (0.2, 0.5),
    "K_I": (5.0, 20.0),
    "gamma_dil": (0.2, 0.5),
    "gamma_lva": (0.7, 1.7),
}


def _d2_predicted_flip_band(condition, tag_variant, t_h, n_draws=2000, seed=0):
    """5th/50th/95th percentile of the D2-predicted ON_only fraction under flip-kinetics
    parameter uncertainty (leak fractions held fixed; see _KINETIC_RANGES note)."""
    rng = np.random.default_rng(seed)
    kf = rng.uniform(*_KINETIC_RANGES["k_flip"], n_draws)
    KI = rng.uniform(*_KINETIC_RANGES["K_I"], n_draws)
    gd = rng.uniform(*_KINETIC_RANGES["gamma_dil"], n_draws)
    gl = rng.uniform(*_KINETIC_RANGES["gamma_lva"], n_draws)
    draws = np.array([
        _d2_predicted_flip(condition, tag_variant, t_h, k_flip=kf[i], K_I=KI[i],
                            gamma_dil=gd[i], gamma_lva=gl[i])
        for i in range(n_draws)
    ])
    return np.percentile(draws, [5, 50, 95])


def summarize(df):
    bad = set(df["pcr_call"].unique()) - VALID_CALLS
    if bad:
        raise ValueError(f"pcr_call has invalid categories: {bad}. Must be one of {VALID_CALLS}")

    gfp_stats = (df.groupby(["tag_variant", "induction", "time_h"])["GFP_OD_normalized"]
                   .agg(["mean", "std", "count"]).reset_index())

    pcr_counts = (df.groupby(["tag_variant", "induction", "time_h", "pcr_call"])
                    .size().reset_index(name="n"))
    pcr_pivot = pcr_counts.pivot_table(index=["tag_variant", "induction", "time_h"],
                                        columns="pcr_call", values="n", fill_value=0).reset_index()
    for c in VALID_CALLS:
        if c not in pcr_pivot.columns:
            pcr_pivot[c] = 0
    pcr_pivot["total"] = pcr_pivot[list(VALID_CALLS)].sum(axis=1)
    for c in VALID_CALLS:
        pcr_pivot[f"frac_{c}"] = pcr_pivot[c] / pcr_pivot["total"]

    merged = gfp_stats.merge(pcr_pivot, on=["tag_variant", "induction", "time_h"])

    merged["d2_predicted_ON_only"] = merged.apply(
        lambda r: _d2_predicted_flip(r["induction"], r["tag_variant"], r["time_h"]), axis=1)
    merged["d2_residual"] = merged["frac_ON_only"] - merged["d2_predicted_ON_only"]

    bands = merged.apply(
        lambda r: _d2_predicted_flip_band(r["induction"], r["tag_variant"], r["time_h"]), axis=1)
    merged["d2_band_p5"] = bands.apply(lambda b: b[0])
    merged["d2_band_p50"] = bands.apply(lambda b: b[1])
    merged["d2_band_p95"] = bands.apply(lambda b: b[2])
    merged["inside_d2_90pct_band"] = (
        (merged["frac_ON_only"] >= merged["d2_band_p5"]) &
        (merged["frac_ON_only"] <= merged["d2_band_p95"]))

    return merged


def make_comparison_plot(df, outpath="m2_pilot_comparison.png"):
    conditions = ["glucose", "no_inducer", "arabinose"]
    variants = sorted(df["tag_variant"].unique())
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))

    # Panel 1: GFP/OD leak-vs-switching bar chart
    ax = axes[0]
    latest_t = df.groupby(["tag_variant", "induction"])["time_h"].max().reset_index()
    plot_df = df.merge(latest_t, on=["tag_variant", "induction", "time_h"])
    width = 0.35
    x = np.arange(len(conditions))
    for i, v in enumerate(variants):
        vals = [plot_df[(plot_df.tag_variant == v) & (plot_df.induction == c)]["mean"].mean()
                for c in conditions]
        ax.bar(x + i * width, vals, width, label=v)
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(conditions)
    ax.set_ylabel("GFP/OD600 (normalized)")
    ax.set_title("M2 GFP output by condition")
    ax.legend()

    # Panel 2: ON-sequence detection frequency (PCR ON_only fraction)
    ax = axes[1]
    for i, v in enumerate(variants):
        vals = [plot_df[(plot_df.tag_variant == v) & (plot_df.induction == c)]["frac_ON_only"].mean()
                for c in conditions]
        ax.bar(x + i * width, vals, width, label=v)
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(conditions)
    ax.set_ylabel("fraction of colonies: ON_only PCR call")
    ax.set_title("M2 ON-sequence detection frequency")
    ax.legend()

    fig.suptitle("M2 pilot: leak (glucose/no-inducer) vs induced switching (arabinose)")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    return fig


def main(csv_path):
    df = pd.read_csv(csv_path, comment="#")
    summary = summarize(df)
    summary.to_csv("m2_pilot_summary.csv", index=False)
    make_comparison_plot(summary)
    print(summary[["tag_variant", "induction", "time_h", "mean", "frac_ON_only", "frac_OFF_only", "frac_mixed"]]
          .to_string(index=False))
    n_covered = summary["inside_d2_90pct_band"].sum()
    n_total = len(summary)
    print(f"\nDry Model D2 coverage: {n_covered}/{n_total} conditions fall inside the"
          f" D2 90% prediction band (flip-kinetics uncertainty only, leak fractions fixed)")
    print("Saved m2_pilot_summary.csv and m2_pilot_comparison.png")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "M2_pilot_template.csv")
