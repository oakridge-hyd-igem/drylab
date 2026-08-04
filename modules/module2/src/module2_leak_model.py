"""
BioChronos Module 2 — leakiness / false-positive model.

Deployment-timescale question (distinct from Module 3's 24 h pulse-capture
question): with NO signal present, what is the cumulative probability that leak
alone flips the irreversible switch over a multi-day/week deployment, and which
leak-suppression strategy keeps that false-positive rate below a target
(default <1%)?

Because the flip is irreversible, the relevant quantity is the CUMULATIVE flip
probability integrated over the whole deployment window, not the instantaneous
integrase level (Module 2 plan). A small leak over weeks still corrupts the
record.

Shared parameters (k_flip, K_I, gamma_dil, gamma_LVA, beta_leak baseline) are
inherited from Module 3 (module3_flip_model.py / module3_parameter_spec.md v6,
Hsiao et al. 2016-anchored) so the two modules stay numerically consistent.

Model (leak only, A = 0, so beta_max / Hill terms are irrelevant):
    dI/dt = beta_leak - gamma * I,          I(0) = 0
        => I(t) = (beta_leak/gamma) * (1 - exp(-gamma*t))
    per-cell flip hazard  lambda(I) = k_flip * I / (I + K_I)   (Hsiao-anchored)
    cumulative flip prob  F(T) = 1 - exp(-∫_0^T lambda(I(t)) dt)

Designs act on the two leak levers:
    weak RBS  -> lower leak translation  (beta_leak scaled by f_rbs)
    LVA tag   -> faster clearance        (gamma += gamma_LVA)
    combined  -> both
"""
import numpy as np

# --- Baselines inherited from Module 3 spec v6 (Hsiao et al. 2016) ---
BASE = dict(
    k_flip=0.4,        # 1/h   max per-cell flip hazard (Hsiao fitted)
    K_I=10.0,          # copies integrase at half-max flip (Hsiao Kd)
    gamma_dil=0.3,     # 1/h   untagged clearance (Hsiao k_deg)
    gamma_LVA=1.04,    # 1/h   extra clearance from LVA ssrA tag (Andersen 1998)
    beta_leak0=100.0,  # copies/cell/h  untagged, strong-RBS leak. GROUNDED (T1, 2026-07-27):
    #                  10% of beta_max, from measured natural CueR/P(copA) copper-promoter
    #                  fold-change 7-18x (Fu 2024) -> leak fraction 5.6-14%, baseline 10%.
    #                  Module 7 sweep 10-500 (1% engineered floor to 50% worst natural).
    #                  [Previous value 10.0 (1% Hsiao integrase construct) was optimistic for copper.]
    f_rbs_weak=0.1,    # weak RBS translation factor (10x weaker; matches Module 3 beta_leak 10->1)
)

DESIGNS = {
    "untagged": dict(f_rbs=1.0,  lva_tag=False),
    "weak RBS": dict(f_rbs=0.1,  lva_tag=False),
    "LVA tag":  dict(f_rbs=1.0,  lva_tag=True),
    "combined": dict(f_rbs=0.1,  lva_tag=True),
}


def design_params(design, base=BASE):
    """(beta_leak, gamma) for a named/dict design under no-signal conditions."""
    d = DESIGNS[design] if isinstance(design, str) else design
    beta_leak = base["beta_leak0"] * d["f_rbs"]
    gamma = base["gamma_dil"] + (base["gamma_LVA"] if d["lva_tag"] else 0.0)
    return beta_leak, gamma


def hazard(I, base=BASE):
    """Per-cell flip hazard, Hsiao-anchored saturating (Michaelis-Menten) form."""
    return base["k_flip"] * I / (I + base["K_I"])


def I_leak(t, beta_leak, gamma):
    """Deterministic leak integrase accumulation, I(0)=0."""
    return (beta_leak / gamma) * (1.0 - np.exp(-gamma * t))


def cumulative_fp(T_grid, beta_leak, gamma, base=BASE, dt=None):
    """Expected leak-driven flipped fraction F(T) = 1 - exp(-∫ lambda dt).

    Analytic mean-field result (validated against the Module 3 agent-based
    stochastic engine to <=0.008 at 24 h). T_grid in hours.
    """
    T_grid = np.atleast_1d(np.asarray(T_grid, float))
    tmax = float(T_grid.max())
    if dt is None:
        dt = min(0.005, tmax / 200000) if tmax > 0 else 0.005
    tt = np.arange(0, tmax + dt, dt)
    lam = hazard(I_leak(tt, beta_leak, gamma), base)
    Lam = np.concatenate([[0.0], np.cumsum(0.5 * (lam[1:] + lam[:-1]) * dt)])
    F = 1.0 - np.exp(-Lam)
    return np.interp(T_grid, tt, F)


def time_to_fp(F_thresh, beta_leak, gamma, base=BASE, t_max_h=24 * 7 * 8):
    """Deployment time (h) at which cumulative FP first crosses F_thresh."""
    tt = np.arange(0, t_max_h, 0.01)
    F = cumulative_fp(tt, beta_leak, gamma, base)
    idx = np.argmax(F >= F_thresh)
    return tt[idx] if F[idx] >= F_thresh else np.inf


def required_leak_integrase(T_hours, F_thresh=0.01, base=BASE):
    """Max steady-state leak integrase (copies/cell) to keep FP < F_thresh over T.

    Plateau approximation Lam(T) ~ lambda(I_ss)*T (transient negligible over
    days). Inverts the saturating hazard for I_ss.
    """
    lam_max = -np.log(1 - F_thresh) / T_hours
    kf, KI = base["k_flip"], base["K_I"]
    if lam_max >= kf:
        return 0.0  # unreachable: even zero-K_I saturating rate exceeds budget
    return lam_max * KI / (kf - lam_max)


def rank_designs(readout_windows=(24.0, 168.0, 672.0), F_thresh=0.01, base=BASE):
    """Summary table across the four designs: leak integrase, hazard, time-to-threshold."""
    rows = []
    base_lam = None
    for name in DESIGNS:
        bl, g = design_params(name, base)
        Iss = bl / g
        lam_ss = hazard(Iss, base)
        if base_lam is None:
            base_lam, base_Iss = lam_ss, Iss
        row = dict(design=name, beta_leak=bl, gamma=g, I_leak_ss=Iss,
                   lambda_ss=lam_ss,
                   rate_fold_reduction=base_lam / lam_ss,
                   Iss_fold_reduction=base_Iss / Iss,
                   t_to_thresh_min=time_to_fp(F_thresh, bl, g, base) * 60,
                   FP_24h=cumulative_fp(24.0, bl, g, base)[0])
        rows.append(row)
    return rows


if __name__ == "__main__":
    for r in rank_designs():
        print(f"{r['design']:9} I_leak_ss={r['I_leak_ss']:7.3f}  "
              f"{r['Iss_fold_reduction']:5.1f}x   "
              f"t(1%FP)={r['t_to_thresh_min']:5.1f} min   FP@24h={r['FP_24h']:.3f}")
