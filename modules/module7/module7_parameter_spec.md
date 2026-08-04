# Module 7 - Sensitivity Analysis (OAT)

**BioChronos (Oakridge-HYD iGEM 2026) | analyte: COPPER | 2026-07-28**

Maps to the plan's Module 7 (per-module sensitivity) and to Nishida-sensei's
"which parameters are worth defending / measuring" question.

## 1. Method

One-at-a-time (OAT) sweep: each parameter is moved from its low to its high plausible
bound while every other parameter stays at its grounded baseline. For each sweep we
record the resulting value of four headline metrics and the **maximum absolute relative
swing** vs the baseline. A large swing means the conclusion depends on that parameter, so
it must be measured precisely and defended; a small swing means the conclusion is robust
to it. Metrics (all for the recommended leak-suppressed **combined** design):

- **A. Trust window** = time for leak-only cumulative false-positive to reach 1% (min) [Module 2]
- **B. Detection threshold** at a 2 h readout (uM Cu, net signal over leak background) [Module 1]
- **C. 90-generation memory retention** (post-event f0 = 0.90) [Module 4]
- **D. 3-state classification accuracy** at 1.25 h readout [Module 5]

Parameters, baselines, plausible ranges and their sources are in the code (PARAMS table)
and `module7_results.csv`. Ranges are the plan's "plausible range" column where given,
else literature bounds (Hsiao, Andersen, Fu, Taniguchi) or explicit design assumptions.

## 2. Results: what each metric is sensitive to

**A. Trust window (baseline ~15 min):**
| Parameter | Swing | Note |
|---|---|---|
| **beta_leak (promoter leak)** | **231%** | dominant: 52 min at 1% leak to 8 min at 50% leak |
| k_flip (Bxb1 rate) | 46% | faster flip shortens the window |
| K_I (flip half-sat) | 38% | higher K_I lengthens it |
| f_rbs (weak RBS) | 38% | tighter RBS lengthens it |
| gamma_dil, gamma_LVA | ~0% | degradation barely moves leak-only FP timing |

**B. Detection threshold (baseline 0.25 uM):**
| Parameter | Swing | Note |
|---|---|---|
| **beta_leak** | **1540%** | dominant: leak floods the readout, blows up the net threshold |
| **n (Hill coefficient)** | **380%** | n=2 sharpens the response, raising the low-Cu net threshold |
| **K (sensor half-max)** | **196%** | threshold tracks K directly |
| k_flip | 108% | affects how fast signal accumulates |
| f_rbs, K_I, gamma | <=32% | secondary |

**C. 90-generation retention (baseline 100%):**
| Parameter | Swing |
|---|---|
| every parameter | **~0%** |

Memory retention is **insensitive to every parameter across its full range** (retention
stays >=99.9%, low 99.89%). Because the record is a heritable DNA inversion and residual leak only
reinforces it, memory stability is not a binding constraint no matter where the parameters
land. This is the strongest robustness result in the project.

**D. 3-state accuracy (baseline 94%):**
| Parameter | Swing | Note |
|---|---|---|
| **CV (expression noise)** | **36%** | dominant: 99.7% at CV=0.1 to 61% at CV=1.0 (task T3) |
| **beta_leak** | **29%** | high leak collapses the states toward each other |
| K (sensor half-max) | 9% | secondary |
| n (Hill) | 6% | secondary |
| others | <=2% | negligible |

## 3. The parameter priority list (answer to Nishida-sensei)

**Tier 1 - measure precisely, be ready to defend (drive multiple conclusions):**
1. **beta_leak / uninduced flip rate** - the single most important parameter in the whole
   project. It is the top driver of the trust window and the detection threshold, and a
   major driver of multi-state separability. Grounded to Fu 2024 copper-promoter data, but
   the authoritative value is the wet-lab EXP1 uninduced-flip measurement.
2. **Bxb1 flip kinetics (k_flip, K_I)** - drive the trust window and threshold; currently a
   Hsiao 2016 proxy, needs pilot wet-lab fitting (the acknowledged gap).
3. **Copper sensor K and Hill n** - dominate the detection threshold; currently inferred
   from operating windows (K) and the MerR single-site mechanism (n~1). The EXP3
   dose-response fit pins both.

**Tier 2 - confirm, but conclusions are more robust:**
4. **Expression-noise CV** - only matters for the multi-state readout (metric D); the
   binary recorder does not depend on it. Task T3 (measure in EXP3).

**Tier 3 - do not need precise values (swept, conclusions robust):**
5. **Fitness cost c, growth/dilution gamma, LVA degradation gamma_LVA, weak-RBS f_rbs** -
   none materially move any headline conclusion across its plausible range. In particular,
   memory retention (metric C) is flat against all of them.

## 4. What this means for the design

- The recorder's performance is **governed by leak** (tier 1, parameter 1) far more than by
  any other quantity. This reinforces every prior module: leak control is the binding
  constraint, and pinning down the real uninduced flip rate is the highest-value wet-lab
  measurement.
- The **sensor parameters (K, n)** matter specifically for the detection threshold, which is
  why the EXP3 dose-response is the second measurement priority.
- **Memory stability is not a risk** under any parameter combination - the team can
  de-prioritise EXP2 with confidence.
- The **multi-state readout** is the only place expression noise matters; if the binary
  recorder is the deliverable, noise is not on the critical path.

## 5. Caveats

1. **OAT, not global.** This is a one-at-a-time analysis; it does not capture interactions
   between parameters (e.g. leak x k_flip jointly). For the design questions here the
   dominant effects are individually large enough that OAT gives the right priorities, but a
   full variance-based (Sobol) analysis would be the rigorous follow-up if two tier-1
   parameters prove strongly coupled.
2. **Metrics are the combined design.** The rankings are for the recommended leak-suppressed
   design; the untagged design is more leak-dominated still.
3. **Ranges are plausibility bounds, not measured error bars.** Once the wet lab returns
   EXP1/EXP3 values, the tier-1 ranges should be replaced with measured confidence intervals
   and the sweep re-run.

## 6. References / sources

Parameter baselines and ranges inherit from Modules 1-5 and their sources: Fu et al. 2024
(*Commun. Biol.* 7:1407, copper-promoter leak); Hsiao et al. 2016 (*Mol. Syst. Biol.* 12:869,
flip kinetics); Andersen et al. 1998 (*Appl. Environ. Microbiol.* 64:2240, ssrA-LVA);
Taniguchi et al. 2010 (*Science* 329:533, expression-noise floor). Metric implementations
re-use the same flip engine (first-order hazard) validated in Module 3.
