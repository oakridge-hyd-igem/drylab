# BioChronos — Module 3: Stochastic Flip Model

**Dry lab / modelling** · iGEM

A biological event recorder. Engineered *E. coli* converts a transient water-contamination
event into a **permanent, irreversible DNA flip** (Bxb1 integrase, attB/attP recombination),
read out later by flow cytometry or a PCR orientation assay. The core claim is **event
recovery** — recording that something *happened*, not measuring the present.

Module 3 is the stochastic model of the flip itself. It answers one question:

> **What is the minimum contamination pulse that leaves a detectable (≥5%) permanent record —
> and what limits it?**

---

## TL;DR — the headline result

**The recorder is *leak-limited*, not *sensitivity-limited*.**

Because the flip is **irreversible** and integrase persists for ~1/γ after a pulse ends, the
device *integrates* flip signal over time. This has one dominant consequence: given enough
time, even tiny basal ("leak") integrase expression flips the whole population **with no
contamination at all** — a false positive. So the practical limit on the device is not how
*short* a pulse it can catch, but how *long* it stays trustworthy before leak trips it.

We define the **trust horizon**: the time before leak alone crosses the 5% detection
threshold. Leak-suppression designs are compared by how much trust horizon they buy.

| Design | γ (1/h) | leak integrase (copies) | **Trust horizon** | False-pos @ 24 h |
|---|---|---|---|---|
| untagged | 0.30 | 33.3 | 0.61 h | 100% |
| LVA tag | 1.34 | 7.5 | 0.67 h | 98% |
| weak RBS | 0.30 | 3.3 | 1.84 h | 88% |
| **combined (weak RBS + LVA)** | 1.34 | 0.75 | **2.59 h** | 47% |

**Findings:**
1. The **combined** design gives the longest trust horizon — **~4× the untagged baseline** —
   by cutting steady-state leak integrase ~44×.
2. Within the trust window, **sub-minute to ~2-minute pulses are recordable** if read at the
   right time. Sensitivity is not the binding constraint; specificity (leak) is.
3. The same leak suppression that fixes false positives costs almost nothing in short-pulse
   sensitivity — the two failure modes are in tension, but at these baselines the tension
   resolves in favour of suppressing leak.

![Module 3 results]({{artifact:art_0f3915f2-eec0-44cc-9767-90dc73d34ee0}})

*Panels: (a) an irreversible flip records a 1 h pulse, but leak alone eventually flips the
population too; (b) leak suppression extends the trust horizon; (c) within the trust window,
short pulses remain recordable if read late enough; (d) trust horizon by design.*

---

## The model

Per-cell state: integrase count `I(t)` and flip status `s ∈ {0,1}` (latches at 1,
irreversible). Under the project QSSA convention (mRNA at steady state) we model integrase
**protein** directly. Production and clearance are deterministic; **only the flip is
stochastic**, because the readout is the population flipped fraction of a discrete,
irreversible per-cell event.

```
Analyte pulse:      A(t) = A_pulse for 0 ≤ t < t_pulse, else 0

Integrase (per cell, deterministic Euler):
    dI/dt = β(A) − γ·I
    β(A)  = β_leak + β_max · Aⁿ / (Kⁿ + Aⁿ)        # Hill activation + basal leak
    γ     = γ_dil                                    # untagged
    γ     = γ_dil + γ_LVA                            # with LVA degradation tag

The flip (stochastic, per unflipped cell — Poisson hazard):
    λ(I) = k_flip · I / (I + K_I)                    # saturating in integrase
    P(flip in [t, t+dt] | s=0) = 1 − exp(−λ(I)·dt)

Readout: event recorded if flipped fraction F = (1/N)Σsᵢ ≥ F_detect
```

---

## Parameters (baselines)

All parameters live in a single `config` dict (`default_config()`), so a one-at-a-time
sensitivity sweep (Module 7) can wrap `simulate()` without refactoring. Baselines are anchored
to the fitted parameters of the one published Bxb1 **recording** model (Hsiao et al. 2016),
with textbook *E. coli* values (Uri Alon) inside the sweep ranges.

| Symbol | Meaning | Baseline | Range (sweep) | Source |
|---|---|---|---|---|
| `β_max` | max integrase production | 1000 copies·cell⁻¹·h⁻¹ | 100–10 000 | Alon, Ch. 2 (est.) |
| `n` | Hill coefficient | 2 | 1–4 | Alon, Ch. 2 |
| `K` | half-max analyte | 1 µM | 0.1–10 | sensor TBD (est.) |
| `β_leak` | basal production (leak) | 10 (=1%·β_max) | 1–100 | Hsiao 2016 (1–2%; measured 0.5–3%) |
| `γ_dil` | untagged clearance | 0.3 h⁻¹ | 0.3–2.1 | Hsiao 2016 (k_deg) |
| `γ_LVA` | extra clearance, LVA tag | 1.04 h⁻¹ | 0.7–1.7 | Andersen et al. 1998 |
| `k_flip` | max flip hazard | 0.4 h⁻¹ | 0.1–10 | Hsiao 2016; Zhang et al. 2022 (tunable) |
| `K_I` | integrase at half-max flip | 10 copies | 1–1000 | Hsiao 2016 (Kd) |
| `A_pulse` | analyte during pulse (input) | 10 µM (~10·K) | 0–100 | input |
| `t_pulse` | pulse duration (input) | scanned | 0.05–24 h | input |
| `F_detect` | detection threshold | 0.05 | 0.01–0.50 | flow cytometry |
| `N` | population size | 2000 | 100–100 000 | numerical |
| `dt` | timestep | 0.01 h | 0.001–0.05 | numerical |

Least-certain parameters (Module 7 priorities): **`k_flip`, `β_leak`, `K_I`** — they sit
directly on the false-positive / false-negative tension.

---

## Controls

| Condition | Setting | F_final | Interpretation |
|---|---|---|---|
| zero-analyte, leak on | `t_pulse=0, β_leak=10` | 1.00 | false positives (leak flips) |
| zero-leak, no pulse | `t_pulse=0, β_leak=0` | **0.00** | true negative |
| zero-leak, 1 h pulse | `t_pulse=1, β_leak=0` | 1.00 | pulse recorded without leak |
| baseline, 1 h pulse | `t_pulse=1, β_leak=10` | 1.00 | true event |

The zero-leak/no-pulse control gives F = 0, confirming the false positives come specifically
from leak, not a simulation artefact.

---

## Assumptions & limitations

- **Extrinsic noise not modelled.** Integrase `I(t)` is a single deterministic trajectory
  shared across cells; all stochasticity is in the flip. Cell-to-cell variation in integrase
  *expression* would broaden the flipped-fraction spread — a planned extension.
- **Flip hazard is simplified.** Hsiao et al. model the flip via a tetramer (4th-order in
  integrase). We use a first-order saturating hazard `λ(I) = k_flip·I/(I+K_I)` — appropriate
  for the model's scope, but less steeply integrase-dependent near threshold.
- **Results are parameter-dependent.** The "leak-limited" conclusion holds *at
  literature-anchored baselines*. Module 7 quantifies how far it moves under the sweep.
- **Sensor not finalised** (Cu / nitrate); `K`, `n` are placeholders. `A_pulse` is set
  saturating (~10·K) so the minimum-pulse result is insensitive to them.

---

## Usage

```bash
python module3_flip_model.py     # runs a baseline self-test
```

```python
from module3_flip_model import default_config, simulate

cfg = default_config()
cfg["t_pulse"] = 1.0             # 1-hour contamination pulse
cfg["lva_tag"] = True            # enable LVA degradation tag
r = simulate(cfg)
print(r["F_final"], r["detected"])
```

---

## Files in this repo

| File | Description |
|---|---|
| `module3_flip_model.py` | The simulation — `default_config()` + `simulate(cfg)`, config-dict convention |
| `module3_parameter_spec.md` | Full parameter specification: equations, sourced values, flagged uncertainties |
| `module3_results.csv` | Per-design trust horizon, leak integrase, false-positive fraction |
| `module3_results.png` | Four-panel results figure |
| `module3_verification.md` | Source verification: β_leak, Hsiao citation, controls |

---

## References

- U. Alon, *An Introduction to Systems Biology*, Ch. 1–2 — baseline *E. coli* production,
  dilution (γ = ln2/τ), Hill parameters.
- V. Hsiao, Y. Hori, P. W. K. Rothemund, R. M. Murray (2016). *A population-based temporal
  logic gate for timing and recording chemical events.* **Mol. Syst. Biol.** 12:869.
  doi:10.15252/msb.20156663
- Q. Zhang, S. M. Azarin, C. A. Sarkar et al. (2022). *Model-guided engineering of DNA
  sequences with predictable site-specific recombination rates.* **Nat. Commun.** 13:4152.
  doi:10.1038/s41467-022-31538-3
- J. B. Andersen, C. Sternberg, L. K. Poulsen et al. (1998). *New unstable variants of green
  fluorescent protein for studies of transient gene expression in bacteria.* **Appl. Environ.
  Microbiol.** 64:2240–2246.
