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
    - m2_pilot_summary.csv    : per (tag_variant, induction, time_h) GFP stats + PCR category fractions
    - m2_pilot_comparison.png : GFP/OD600 leak and switching, untagged vs LVA, by condition
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

VALID_CALLS = {"OFF_only", "mixed", "ON_only"}


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
    print("Saved m2_pilot_summary.csv and m2_pilot_comparison.png")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "M2_pilot_template.csv")
