"""
BioChronos Module 1 — Copper sensor dose-response & recording threshold.

Analyte: COPPER (locked 2026-07-27). Sensor: CueR / P(copA) — the Cue system,
the primary and best-characterised E. coli copper sensor (MerR-family activator,
aerobic, cytoplasmic Cu+).

Delivers (per the modelling plan's Module 1 / Sensor Threshold):
  - integrase production dose-response  beta(Cu)
  - recording dose-response  F(Cu)  (flip fraction vs copper), coupled through
    the Module 3 flip engine so Models 1/2/3 stay numerically consistent
  - activation / detection threshold (min recordable copper)
  - Hill coefficient  (grounded: n ~= 1, near-hyperbolic — NOT cooperative)

Literature-grounded parameters (see module1_parameter_spec.md for full sourcing):
  n (Hill)      ~ 1.0    Stoyanov 2001 abstract (copA::lacZ "approximately
                         proportional to the concentration of cupric ions")
                         + MerR-family single-site mechanism -> non-cooperative.
                         (An earlier draft cited a Zahid 2012 "Hill=1.0" figure;
                         that figure is WITHDRAWN — full text not retrievable.)
  K (half-max)  ~ 5 uM   inferred from operating windows (no explicit fitted K
                         in the literature); grounded plausible range 1-15 uM
  LOD           0.15-0.25 uM   Pang 2020 (WMC-007, 0.25); Ivask 2009 (0.15)
  linear range  ~0.4-80 uM     Pang 2020 (0.39-78.68); Riether 2001 (3-30)
  regulatory    WHO 30 uM (2 mg/L), China 15 uM (1 mg/L)   Fu 2024

Flip kinetics (k_flip, K_I) and leak/clearance are inherited unchanged from
Module 3 (Hsiao 2016-anchored).
"""
import numpy as np
import importlib.util

_trapz = getattr(np, "trapezoid", None) or np.trapz  # numpy>=2 renamed trapz->trapezoid

# --- grounded copper-sensor baselines ---
K_COPPER   = 5.0    # uM, half-max (plausible range 1-15)
N_HILL     = 1.0    # near-hyperbolic (grounded)
BETA_MAX   = 1000.0 # copies/cell/h, fully-induced production (Module 3)
BETA_LEAK  = 100.0  # copies/cell/h, basal. GROUNDED (T1, 2026-07-27): 10% of beta_max,
#                    from measured natural CueR/P(copA) copper-promoter fold-change 7-18x
#                    (Fu 2024) -> leak fraction 5.6-14%, baseline 10%. Was 10.0 (1% Hsiao,
#                    optimistic for copper). Module 7 sweep 10-500 (1%-50%).
K_FLIP     = 0.4    # 1/h  (Module 3 / Hsiao)
K_I        = 10.0   # copies (Module 3 / Hsiao)
GAMMA_DIL  = 0.3    # 1/h  (Module 3)
GAMMA_LVA  = 1.04   # 1/h  (Module 3)
F_DETECT   = 0.05   # 5% flow-cytometry detectable fraction (Module 3)

# Deployment designs (from Module 2): untagged baseline vs combined leak-suppressed
DESIGNS = {
    "untagged": dict(gamma=GAMMA_DIL,            beta_max=BETA_MAX,      beta_leak=BETA_LEAK),
    "combined": dict(gamma=GAMMA_DIL+GAMMA_LVA,  beta_max=BETA_MAX*0.1,  beta_leak=BETA_LEAK*0.1),  # weak RBS x0.1
}


def beta_of_Cu(Cu, beta_max=BETA_MAX, beta_leak=BETA_LEAK, K=K_COPPER, n=N_HILL):
    """Integrase production rate (copies/cell/h) as a Hill function of copper (uM)."""
    Cu = np.asarray(Cu, float)
    return beta_leak + beta_max * Cu**n / (K**n + Cu**n)


def hazard(I, k_flip=K_FLIP, K_I=K_I):
    """Per-cell flip hazard (Hsiao-anchored saturating form)."""
    return k_flip * I / (I + K_I)


def flip_fraction(Cu, gamma, t_expose, k_flip=K_FLIP, K_I=K_I, **kw):
    """Flip fraction after sustained copper exposure of length t_expose (h).

    Integrase ramps from 0 to steady state during exposure; cumulative flip
    probability F = 1 - exp(-int lambda dt). Vectorised over Cu.
    """
    Cu = np.atleast_1d(np.asarray(Cu, float))
    tt = np.arange(0, t_expose + 1e-9, 0.01)
    Iss = beta_of_Cu(Cu, **kw) / gamma
    I = Iss[:, None] * (1 - np.exp(-gamma * tt)[None, :])
    Lam = _trapz(hazard(I, k_flip, K_I), tt, axis=1)
    return 1 - np.exp(-Lam)


def detection_threshold(design, t_expose, margin=F_DETECT, K=K_COPPER, n=N_HILL,
                        Cu_grid=None):
    """Net signal-over-background copper threshold (uM): min Cu adding >=margin
    flip fraction ABOVE the no-signal leak background. This is the correct
    activation threshold for an irreversible recorder with non-zero leak."""
    if Cu_grid is None:
        Cu_grid = np.logspace(-2, 3, 4000)
    d = DESIGNS[design]
    F = flip_fraction(Cu_grid, d["gamma"], t_expose, K=K, n=n,
                      beta_max=d["beta_max"], beta_leak=d["beta_leak"])
    F0 = flip_fraction([0.0], d["gamma"], t_expose, K=K, n=n,
                       beta_max=d["beta_max"], beta_leak=d["beta_leak"])[0]
    net = F - F0
    idx = np.argmax(net >= margin)
    thr = Cu_grid[idx] if net[idx] >= margin else np.inf
    return thr, F0


if __name__ == "__main__":
    print("Net copper detection threshold (uM) to add >=5% flips above leak:")
    for te in [1, 2, 3, 4]:
        tu, bu = detection_threshold("untagged", te)
        tc, bc = detection_threshold("combined", te)
        print(f"  {te}h  untagged={tu:6.3f} (leak bg {bu*100:4.0f}%)   "
              f"combined={tc:6.3f} (leak bg {bc*100:4.0f}%)")
