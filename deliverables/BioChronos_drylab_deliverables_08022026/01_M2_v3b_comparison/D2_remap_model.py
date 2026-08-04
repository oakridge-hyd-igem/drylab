"""
BioChronos Dry Model D2 (re-mapped) - weak-RBS Untagged vs weak-RBS+LVA
=======================================================================
Section 2 re-mapping of the M2 v3b Revision Instructions.

Isolated-variable comparison: the two variants differ ONLY in the Bxb1
degradation term (LVA tag -> +gamma_LVA). RBS, copy number, promoter, and
flip kinetics are identical. This is the weak-RBS-untagged vs same-weak-RBS
+ LVA comparison the wet lab is actually building (Wet EXP-2), NOT the old
strong-RBS-untagged vs combined pairing.

Leak is grounded to the pBAD/AraC arabinose promoter (current M2 is a pBAD
memory test), NOT the copper promoter. Copper-promoter leak (Fu 2024,
beta_leak~100) is reserved for the FUTURE integrated copper recorder (Wet M3).

pBAD leak fractions (promoter-set, RBS-independent):
  glucose repression : 1/1200  (Guzman et al. 1995, J Bacteriol 177:4121-4130;
                                induced/repressed ratio up to 1200-fold)     [literature proxy]
  no inducer         : 0.005    (leakier OFF; no arabinose, no glucose;
                                bracket 0.1-2%, baseline 0.5%)               [assumed]
  arabinose (ON)     : 1.0      (defines beta_max_eff; true-positive)         [defines scale]

Kinetics shared with the grounded engine: k_flip=0.4/h, K_I=10 copies
(Hsiao 2016); gamma_dil=0.3/h; gamma_LVA=1.04/h added for the LVA variant
(Andersen 1998 ssrA-LVA t1/2 ~40 min).
"""
import numpy as np

BETA_MAX_STRONG = 1000.0   # strong-RBS fully-induced integrase production (copies/cell/h)
F_RBS_WEAK      = 0.1      # B0033 weak RBS
BETA_MAX_EFF    = BETA_MAX_STRONG * F_RBS_WEAK   # = 100, induced weak-RBS
K_FLIP, K_I     = 0.4, 10.0
GAMMA_DIL       = 0.30
GAMMA_LVA       = 1.04

LEAK_FRAC = {"glucose repression": 1/1200, "no inducer": 0.005, "arabinose (ON)": 1.0}

def gamma(variant):
    """variant in {'Untagged','LVA'}; LVA adds active degradation."""
    return GAMMA_DIL + (GAMMA_LVA if variant == "LVA" else 0.0)

def beta_of(cond):
    """Production rate under a condition: induced beta_max for ON, leak fraction*beta_max for OFF."""
    return BETA_MAX_EFF if cond == "arabinose (ON)" else LEAK_FRAC[cond] * BETA_MAX_EFF

def cumulative_flip(cond, variant, tmax, dt=1e-3):
    """Flipped fraction over [0,tmax]. I(t) deterministic; per-cell hazard k_flip*I/(I+K_I)."""
    bl, g = beta_of(cond), gamma(variant)
    tt = np.arange(0, tmax + dt, dt)
    I  = (bl / g) * (1 - np.exp(-g * tt))
    lam = K_FLIP * I / (I + K_I)
    cum = np.concatenate([[0], np.cumsum((lam[1:] + lam[:-1]) / 2 * np.diff(tt))])
    return tt, 1 - np.exp(-cum)

def flip_at(cond, variant, t):
    _, F = cumulative_flip(cond, variant, t)
    return F[-1]

if __name__ == "__main__":
    print("beta_max_eff (weak RBS induced) =", BETA_MAX_EFF)
    for cond in LEAK_FRAC:
        for v in ["Untagged", "LVA"]:
            print(f"  {cond:20s} {v:9s} flip@8h={flip_at(cond,v,8)*100:5.1f}%  flip@16h={flip_at(cond,v,16)*100:5.1f}%")
