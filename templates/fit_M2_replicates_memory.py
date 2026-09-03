"""
fit_M2_replicates_memory.py

Analyzes M2 biological replicates and the inducer-removal memory test.
Computes:
  1. Population-to-population CV of GFP/OD (resolves task T3 - decision hinge for
     single- vs two-integrase architecture, per the consolidated report).
  2. Memory decay: GFP/OD and ON-sequence PCR fraction over hours since inducer removal.

Usage:
    python fit_M2_replicates_memory.py M2_replicates_memory_data.csv

Input CSV columns (see M2_replicates_memory_template.csv):
    construct, tag_variant, induction, bio_replicate_n, time_h,
    hours_since_inducer_removal, OD600, GFP_raw, GFP_OD_normalized, pcr_call,
    operator, date

Output:
    - m2_replicate_cv.csv       : CV of GFP/OD per (tag_variant, induction, time_h)
    - m2_memory_decay.csv       : mean GFP/OD and ON_only fraction vs hours since removal
    - m2_replicate_memory.png   : CV summary + memory decay curve
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

CV_DECISION_THRESHOLD = 0.5  # from the consolidated report: CV > ~0.5 favors two-integrase


def compute_cv(df):
    grp = df.groupby(["tag_variant", "induction", "time_h"])["GFP_OD_normalized"]
    cv = grp.agg(lambda s: s.std(ddof=1) / s.mean() if s.mean() != 0 else np.nan)
    n = grp.count()
    out = pd.DataFrame({"CV": cv, "n_bio_replicates": n}).reset_index()
    out["exceeds_0.5_threshold"] = out["CV"] > CV_DECISION_THRESHOLD
    return out


def compute_memory_decay(df):
    mem = df.dropna(subset=["hours_since_inducer_removal"])
    if mem.empty:
        return pd.DataFrame()
    gfp = (mem.groupby(["tag_variant", "hours_since_inducer_removal"])["GFP_OD_normalized"]
              .mean().reset_index())
    mem = mem.assign(_is_on=(mem["pcr_call"] == "ON_only"))
    on_frac = (mem.groupby(["tag_variant", "hours_since_inducer_removal"])["_is_on"]
                  .mean().reset_index(name="frac_ON_only"))
    return gfp.merge(on_frac, on=["tag_variant", "hours_since_inducer_removal"])


def make_plot(cv_df, decay_df, outpath="m2_replicate_memory.png"):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))

    ax = axes[0]
    labels = cv_df["tag_variant"] + " / " + cv_df["induction"]
    colors = ["#c0392b" if e else "#1f6fb4" for e in cv_df["exceeds_0.5_threshold"]]
    ax.barh(labels, cv_df["CV"], color=colors)
    ax.axvline(CV_DECISION_THRESHOLD, linestyle="--", color="gray",
               label=f"decision threshold ({CV_DECISION_THRESHOLD})")
    ax.set_xlabel("CV of GFP/OD across biological replicates")
    ax.set_title("Population-to-population noise (task T3)")
    ax.legend(fontsize=8)

    ax = axes[1]
    if not decay_df.empty:
        for v in decay_df["tag_variant"].unique():
            sub = decay_df[decay_df.tag_variant == v].sort_values("hours_since_inducer_removal")
            ax.plot(sub["hours_since_inducer_removal"], sub["frac_ON_only"], "-o", label=v)
        ax.set_xlabel("hours since inducer removal")
        ax.set_ylabel("fraction of colonies: ON_only PCR call")
        ax.set_title("Memory retention after inducer removal")
        ax.legend()
    else:
        ax.text(0.5, 0.5, "no memory-test rows in data", ha="center")

    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    return fig


def main(csv_path):
    df = pd.read_csv(csv_path, comment="#")
    cv_df = compute_cv(df)
    cv_df.to_csv("m2_replicate_cv.csv", index=False)
    decay_df = compute_memory_decay(df)
    if not decay_df.empty:
        decay_df.to_csv("m2_memory_decay.csv", index=False)
    make_plot(cv_df, decay_df)
    print(cv_df.to_string(index=False))
    if (cv_df["exceeds_0.5_threshold"]).any():
        print("\nCV exceeds 0.5 in at least one condition: two-integrase architecture is favored (per report decision rule).")
    else:
        print("\nCV stays below 0.5: single leak-suppressed integrase remains sufficient (per report decision rule).")
    print("Saved m2_replicate_cv.csv, m2_memory_decay.csv (if applicable), and m2_replicate_memory.png")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "M2_replicates_memory_template.csv")
