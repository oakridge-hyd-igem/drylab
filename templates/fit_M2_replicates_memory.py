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
    - m2_memory_decay.csv       : mean GFP/OD and ON_only fraction vs hours since removal,
                                  plus the Module 4 retention-model prediction band
                                  (fitness cost swept over its grounded 0-5% range)
    - m2_replicate_memory.png   : CV summary + memory decay curve with model band
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

CV_DECISION_THRESHOLD = 0.5  # from the consolidated report: CV > ~0.5 favors two-integrase

# --- Module 4 retention model, reproduced inline (module4_memory_model.py) ---
# df/dt = (1-f)*(lam - MU*c*f); lam = leak hazard AFTER inducer removal (pBAD no-inducer
# leak, not the arabinose-ON rate), c = fitness cost of the flipped state, swept over
# Module 4's grounded range [0, 0.05] (Bxb1 burden + Bienick 2014 tradeoff bound).
_MU = 0.3          # 1/h growth = dilution rate
_K_FLIP, _K_I = 0.4, 10.0
_BETA_MAX_EFF = 1000.0 * 0.1
_NO_INDUCER_LEAK_FRAC = 0.005
_GAMMA_DIL, _GAMMA_LVA = 0.30, 1.04
_COST_RANGE = (0.0, 0.05)


def _leak_hazard_after_removal(tag_variant):
    gamma = _GAMMA_DIL + (_GAMMA_LVA if tag_variant == "LVA" else 0.0)
    beta_leak = _NO_INDUCER_LEAK_FRAC * _BETA_MAX_EFF
    Iss = beta_leak / gamma
    return _K_FLIP * Iss / (Iss + _K_I)


def _retention_curve(t_h, f0, lam, c):
    from scipy.integrate import odeint
    def dfdt(f, t):
        return (1 - f) * (lam - _MU * c * f)
    t_h = np.atleast_1d(t_h).astype(float)
    grid = np.unique(np.concatenate([[0.0], t_h]))
    sol = odeint(dfdt, f0, grid).ravel()
    return np.interp(t_h, grid, sol)


def _module4_prediction_band(t_hours, f0, tag_variant, n_draws=200, seed=0):
    """5th/50th/95th percentile retention curve over the grounded cost range."""
    rng = np.random.default_rng(seed)
    lam = _leak_hazard_after_removal(tag_variant)
    cs = rng.uniform(*_COST_RANGE, n_draws)
    curves = np.array([_retention_curve(t_hours, f0, lam, c) for c in cs])
    return np.percentile(curves, [5, 50, 95], axis=0)


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
    merged = gfp.merge(on_frac, on=["tag_variant", "hours_since_inducer_removal"])
    merged = merged.sort_values(["tag_variant", "hours_since_inducer_removal"]).reset_index(drop=True)

    # Module 4 prediction band: f0 = the measured ON_only fraction at t=0 (removal time)
    # for that tag_variant, so the model is anchored to the observed starting point and
    # only the POST-removal decay/reinforcement shape is being tested, not the induction
    # level itself.
    p5_all, p50_all, p95_all = [], [], []
    for v in merged["tag_variant"].unique():
        sub = merged[merged.tag_variant == v]
        t0_rows = sub[sub["hours_since_inducer_removal"] == sub["hours_since_inducer_removal"].min()]
        f0 = t0_rows["frac_ON_only"].iloc[0] if not t0_rows.empty else 0.9
        t_hours = sub["hours_since_inducer_removal"].values
        p5, p50, p95 = _module4_prediction_band(t_hours, f0, v)
        p5_all.extend(p5); p50_all.extend(p50); p95_all.extend(p95)
    merged["model_p5"] = p5_all
    merged["model_p50"] = p50_all
    merged["model_p95"] = p95_all
    merged["inside_module4_90pct_band"] = (
        (merged["frac_ON_only"] >= merged["model_p5"]) &
        (merged["frac_ON_only"] <= merged["model_p95"]))
    return merged


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
        colors = {"untagged": "#ff7f0e", "LVA": "#1f77b4"}
        for v in decay_df["tag_variant"].unique():
            sub = decay_df[decay_df.tag_variant == v].sort_values("hours_since_inducer_removal")
            c = colors.get(v, None)
            ax.plot(sub["hours_since_inducer_removal"], sub["frac_ON_only"], "o",
                     color=c, label=f"{v} (measured)")
            if "model_p50" in sub.columns:
                ax.plot(sub["hours_since_inducer_removal"], sub["model_p50"], "-", color=c, alpha=0.7)
                ax.fill_between(sub["hours_since_inducer_removal"], sub["model_p5"], sub["model_p95"],
                                 color=c, alpha=0.15,
                                 label=f"{v} Module 4 band (cost 0-5%)" if v == list(decay_df["tag_variant"].unique())[0] else None)
        ax.set_xlabel("hours since inducer removal")
        ax.set_ylabel("fraction of colonies: ON_only PCR call")
        ax.set_title("Memory retention: measured vs Module 4 prediction")
        ax.legend(fontsize=7)
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
        n_covered = decay_df["inside_module4_90pct_band"].sum()
        n_total = len(decay_df)
        print(f"Module 4 retention-model coverage: {n_covered}/{n_total} timepoints fall inside"
              f" the model's 90% band (fitness cost swept 0-5%, f0 anchored to measured t=0)\n")
    make_plot(cv_df, decay_df)
    print(cv_df.to_string(index=False))
    if (cv_df["exceeds_0.5_threshold"]).any():
        print("\nCV exceeds 0.5 in at least one condition: two-integrase architecture is favored (per report decision rule).")
    else:
        print("\nCV stays below 0.5: single leak-suppressed integrase remains sufficient (per report decision rule).")
    print("Saved m2_replicate_cv.csv, m2_memory_decay.csv (if applicable), and m2_replicate_memory.png")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "M2_replicates_memory_template.csv")
