"""
BioChronos Module 3 — stochastic flip model (agent-based).

Deterministic per-cell integrase trajectory (QSSA: mRNA at steady state) +
stochastic, irreversible flip. All parameters live in one config dict so the
Module 7 one-at-a-time sensitivity sweep can wrap simulate() without refactoring.

Baselines: module3_parameter_spec.md v5 (Hsiao et al. 2016-anchored).
Sources per parameter are in that spec; every value is sourced there.
"""
import numpy as np


def default_config():
    """Single source of truth for all parameters (spec v5 baselines)."""
    return dict(
        beta_max=1000.0,   # copies/cell/h  max integrase production (Alon)
        n=2.0,             # Hill coefficient (Alon)
        K=1.0,             # uM  half-max analyte (ESTIMATED, sensor TBD)
        beta_leak=10.0,    # copies/cell/h  basal at A=0 (Hsiao ~1% of prod)
        gamma_dil=0.3,     # 1/h  untagged clearance (Hsiao k_deg)
        gamma_LVA=1.04,    # 1/h  extra from LVA tag (Andersen 1998)
        k_flip=0.4,        # 1/h  max per-cell flip hazard (Hsiao)
        K_I=10.0,          # copies  integrase at half-max flip (Hsiao Kd)
        A_pulse=10.0,      # uM  analyte during pulse (~10*K, saturating)
        t_pulse=1.0,       # h   pulse duration (the quantity we solve for)
        t_total=24.0,      # h   total simulated time
        lva_tag=False,     # LVA degradation tag present?
        F_detect=0.05,     # detectable flipped fraction (flow cytometry)
        N=2000,            # population size
        dt=0.01,           # h  timestep (must be << 1/gamma and << 1/k_flip)
        seed=0,
    )


def beta_of_A(A, cfg):
    """QSSA-collapsed integrase production: Hill activation + basal leak."""
    return cfg["beta_leak"] + cfg["beta_max"] * A**cfg["n"] / (cfg["K"]**cfg["n"] + A**cfg["n"])


def gamma_of(cfg):
    """Total first-order integrase clearance for the chosen design."""
    return cfg["gamma_dil"] + (cfg["gamma_LVA"] if cfg["lva_tag"] else 0.0)


def simulate(cfg):
    """Agent-based population simulation.

    Integrase I(t) is deterministic (Euler ODE) and shared across cells; only
    the flip is stochastic and per-cell (Poisson hazard, irreversible latch).
    Returns dict: t, I, F (time series), F_final, detected, cfg.
    """
    rng = np.random.default_rng(cfg["seed"])
    dt, N = cfg["dt"], cfg["N"]
    nsteps = int(round(cfg["t_total"] / dt))
    gamma = gamma_of(cfg)

    t = np.zeros(nsteps + 1)
    I = np.zeros(nsteps + 1)
    F = np.zeros(nsteps + 1)
    flipped = np.zeros(N, dtype=bool)

    I_cur = 0.0
    for k in range(nsteps):
        tk = k * dt
        A = cfg["A_pulse"] if tk < cfg["t_pulse"] else 0.0
        I_cur = max(I_cur + dt * (beta_of_A(A, cfg) - gamma * I_cur), 0.0)
        lam = cfg["k_flip"] * I_cur / (I_cur + cfg["K_I"])   # flip hazard
        p = 1.0 - np.exp(-lam * dt)
        un = ~flipped
        newly = rng.random(un.sum()) < p
        flipped[np.where(un)[0][newly]] = True
        t[k + 1] = tk + dt; I[k + 1] = I_cur; F[k + 1] = flipped.mean()
    return dict(t=t, I=I, F=F, F_final=F[-1], cfg=cfg,
                detected=F[-1] >= cfg["F_detect"])


if __name__ == "__main__":
    r = simulate(default_config())
    print("baseline untagged, 1h pulse: F_final =", round(r["F_final"], 3),
          "detected =", r["detected"])
