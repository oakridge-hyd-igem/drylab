"""
BioChronos Module 4 — Irreversibility / Memory Stability.

Tests the core memory-persistence claim: once the Bxb1 flip has written a state,
that state must survive after the signal disappears — through continued growth,
cell division, and dilution over a readout delay.

Key insight: the flip is DNA-encoded, so it is HERITABLE. Daughters inherit the
flipped chromosome. Unlike a protein/metabolite memory (which halves every
generation once production stops), the DNA record does not dilute with growth.
The only things that can erode it are (i) a fitness cost of the flipped state
(flipped cells outgrown by unflipped), and (ii) continued leak (which actually
REINFORCES the flipped fraction, since the flip is one-way).

Two-subpopulation model (U = unflipped, F = flipped), growth rate mu = dilution
rate gamma_dil from Modules 2/3 for consistency:
    dU/dt = mu*U - lam*U
    dF/dt = mu*(1-c)*F + lam*U
with lam = leak flip hazard (Module 2), c = fitness cost of the flipped state.

The flipped fraction f = F/(U+F) obeys the closed form (validated to machine
precision against the full ODE):
    df/dt = (1-f)*(lam - mu*c*f)
Equilibrium:  f* = min(1, lam/(mu*c))   (f* -> 1 as c -> 0).
Cost-only decay (lam=0): odds f/(1-f) decay as 2^(-c*g) over g generations, so
the memory odds half-life is 1/c generations.
"""
import numpy as np
from scipy.integrate import odeint

MU = 0.3            # 1/h growth = dilution rate (Modules 2/3)
TAU = np.log(2)/MU  # generation time (h)
K_FLIP, K_I = 0.4, 10.0

DESIGNS = {   # leak (beta_leak, gamma) from Module 2. beta_leak0 GROUNDED to 100 (T1,
              # 2026-07-27): 10% of beta_max, natural CueR/P(copA) fold-change 7-18x (Fu
              # 2024). Was 10/1 (1% Hsiao). Higher leak only REINFORCES the flipped record.
    "untagged": dict(beta_leak=100.0, gamma=0.30),
    "weak RBS": dict(beta_leak=10.0,  gamma=0.30),
    "LVA tag":  dict(beta_leak=100.0, gamma=1.34),
    "combined": dict(beta_leak=10.0,  gamma=1.34),
}

def leak_hazard(beta_leak, gamma):
    Iss = beta_leak/gamma
    return K_FLIP*Iss/(Iss+K_I)

def dfdt(f, t, lam, c):
    return (1-f)*(lam - MU*c*f)

def retention(t_h, f0, lam, c):
    """Flipped fraction at time(s) t_h (h) from initial f0 under leak lam, cost c.

    Integrates from t=0; returns the fraction(s) at the requested time(s).
    """
    t_h = np.atleast_1d(t_h).astype(float)
    grid = np.unique(np.concatenate([[0.0], t_h]))
    sol = odeint(dfdt, f0, grid, args=(lam, c)).ravel()
    return np.interp(t_h, grid, sol)

def f_star(lam, c):
    """Long-run equilibrium flipped fraction."""
    return 1.0 if c <= 0 else min(1.0, lam/(MU*c))

def odds_half_life_gen(c):
    """Memory odds half-life in generations (cost-only, no leak)."""
    return np.inf if c <= 0 else 1.0/c

if __name__ == "__main__":
    for nm, d in DESIGNS.items():
        lam = leak_hazard(**d)
        f48 = retention(48.0, 0.90, lam, 0.05)[0]
        print(f"{nm:9} lam={lam:.4f}/h  f(48h|c=5%,f0=0.9)={f48:.3f}  f*(c=5%)={f_star(lam,0.05):.3f}")
    for c in (0.02, 0.05, 0.10):
        print(f"cost {c:.0%}: odds half-life {odds_half_life_gen(c):.0f} gen = {TAU/c:.0f} h")


# ============================================================================
# Stochastic DNA-record model (Wright-Fisher) — added after wet-lab-style review.
# The flip state is an "allele" on the DNA, inherited by daughters. Erosion comes
# only from (i) fitness cost c (flipped cells outgrown) and (ii) finite-population
# genetic drift. Leak is ONE-WAY (irreversible flip), so it reinforces, never erases.
# Parameters grounded in measured data (see module4_grounded_anchors.json):
#   - fitness cost c: baseline 2%, plausible range 0-5% (Bienick 2014 expression-
#     growth tradeoff; Bxb1 burden "not measurable", Gonzalez-Colell 2025).
#   - validation target: memory maintained >=90 generations (Siuti 2013, measured).
# ============================================================================
def wright_fisher(N, f0, c, lam_h, n_gen, reps, seed=0):
    """Stochastic flipped-fraction trajectories over n_gen generations.

    N     founding/effective population size
    f0    initial flipped fraction (from a recording event)
    c     fitness cost of the flipped state (0-0.05 grounded)
    lam_h leak flip hazard (per h); converted to per-generation prob (one-way)
    Returns array (reps, n_gen+1) of flipped fractions.
    """
    r = np.random.default_rng(seed)
    m = 1 - np.exp(-lam_h*TAU)                       # leak flip prob per generation
    F = r.binomial(N, f0, size=reps).astype(float)
    traj = np.zeros((reps, n_gen+1)); traj[:, 0] = F/N
    for g in range(n_gen):
        f = F/N
        p = f*(1-c)/(f*(1-c) + (1-f))                # selection (flipped cost c)
        p = p + (1-p)*m                              # one-way leak U->F
        F = r.binomial(N, np.clip(p, 0, 1))          # drift: resample N offspring
        traj[:, g+1] = F/N
    return traj

if __name__ == "__main__":
    # Validation against Siuti 2013 measured >=90-generation retention:
    for nm in ("untagged", "combined"):
        tr = wright_fisher(10000, 0.9, 0.02, leak_hazard(**DESIGNS[nm]), 90, 500, seed=11)
        print(f"[WF] {nm:9} c=2%: flipped@90gen = {tr[:,-1].mean():.3f}  (Siuti: memory maintained >=90 gen)")
