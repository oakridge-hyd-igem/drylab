# Module 3 — Pre-Module-7 verification memo

Four checks requested before the sensitivity sweep. All resolved; the headline result stands.

---

## 1. β_leak verified against source (highest priority)

**Claim in spec/model:** baseline `β_leak = 10 copies·cell⁻¹·h⁻¹ = 1% of β_max`, attributed to
"Hsiao et al. 2016 fit leak = 1–2% of production."

**What the Hsiao PDF actually says** (extracted from the fetched full text,
10.15252/msb.20156663):

- **Model parameter values:** the fitted/used leak rates are `k_leakA = 1%` and
  `k_leakB = 2%` of the corresponding production rate `k_prod` — verbatim:
  *"k_leakA = 1%, k_leakB = 2%"* and *"k_leakB = 0.02·k_prodB"*. This is the source of the
  "1–2% of production" claim. **Confirmed accurate.**
- **Independent experimental measurement:** leaky expression was directly measured and
  *"quite low, not exceeding ~0.5–3%"* of induced production — so the 1–2% model values sit
  inside the measured band.

**Verdict:** β_leak baseline is correctly sourced. Our 1% choice = `k_leakA`, the lower of
Hsiao's two fitted values and near the middle of the measured 0.5–3% range. The headline
result (leak-limited recorder, trust horizons 0.6–2.6 h) rests on a value that is grounded in
both a fitted model parameter and a direct measurement.

**Caveat for the sweep:** because the measured spread is 0.5–3% and construct/context
dependence is large, `β_leak` remains the parameter to sweep hardest — the spec range
(1–100 copies = 0.1–10% of β_max) brackets both the measured band and worse-case leak.

---

## 2. Hsiao citation — confirmed real and says what's claimed

Confirmed directly from the fetched PDF header/footer **and** CrossRef (matching fields):

| Field | Value (from PDF + CrossRef) |
|---|---|
| Title | *A population-based temporal logic gate for timing and recording chemical events* |
| Authors | Victoria Hsiao, Yutaka Hori, Paul W. K. Rothemund, Richard M. Murray |
| Journal | *Molecular Systems Biology* **12**: 869 (2016) — the PDF footer prints "Molecular Systems Biology 12: 869 \| 2016", so the article number **869** is legitimate (I had earlier over-cautiously dropped it) |
| DOI | 10.15252/msb.20156663 |

**Does it say what we cite it for?** Yes, each claim checks out against the text:
- *"we compare predictions of **Markov model** simulations with laboratory measurements of
  final population distributions"* — supports our "Markov flip model" attribution.
- *"the timing between inputs, and the duration of input **pulses**"* — supports pulse-input use.
- Fitted constants used as our baselines are present: `k_flip = 0.4 h⁻¹`, `k_deg = 0.3 h⁻¹`,
  `Kd = 10 molecules`, leak `= 1–2%`.

**One nuance worth recording** (does not change our baselines): Hsiao models the flip hazard
as a **tetramer** binding form — *"once both attB and attP sites are occupied, form a tetramer
(dimer of dimers) that digests, flips, and re-ligates the DNA"* — i.e. a 4th-order function of
integrase (their Eq. 2), with binding stated as **non-cooperative** (*"no cooperativity has
been observed in Bxb1 or TP901-1 DNA binding"*). Our model uses a simpler first-order
saturating hazard `λ(I) = k_flip·I/(I+K_I)`. This is a deliberate simplification (appropriate
for a high-school-team model), but it means our flip is *less* steeply integrase-dependent
than Hsiao's near threshold. Flagged as a structural assumption for the sweep / write-up.

---

## 3. Zero-leak control — now run separately from the zero-analyte case

**Finding:** in the original Module 3 run, the false-positive test set `t_pulse = 0` with leak
still present (`β_leak = 10`). A **zero-leak** control (`β_leak = 0`) was **not** run
separately. These are different tests and both are needed to attribute the false positives to
leak specifically. Now run:

| Control condition | Setting | F_final | Expected |
|---|---|---|---|
| (i) zero-analyte, leak ON | `t_pulse=0, β_leak=10` | **1.000** | false positives (leak flips) |
| (ii) **zero-LEAK, no pulse** | `t_pulse=0, β_leak=0` | **0.000** | true negative (nothing flips) |
| (iii) zero-LEAK, 1 h pulse | `t_pulse=1, β_leak=0` | **0.999** | pulse alone still records |
| (iv) baseline, 1 h pulse | `t_pulse=1, β_leak=10` | **1.000** | true event |

**Verdict:** the two tests are distinct and give distinct results. The zero-leak/no-pulse
control gives **F = 0**, confirming the false positives in condition (i) come **specifically
from leak**, not from a simulation artefact. Condition (iii) confirms the flip machinery
records a pulse even with leak fully removed. Control code lives in the analysis; the switch is
already in the config (`β_leak`), so Module 7 can sweep it without changes.

---

## 4. Both caveats retained for the write-up

1. **Extrinsic noise not modelled.** Integrase `I(t)` is a single deterministic trajectory
   shared across cells; all stochasticity is in the flip. Cell-to-cell variation in integrase
   *expression* (extrinsic noise) would broaden the flipped-fraction spread and slightly lower
   the effective threshold-crossing time. Natural Module 7 / extension knob.

2. **Result is parameter-dependent.** The "leak-limited, not sensitivity-limited" conclusion
   holds at the Hsiao-anchored baselines. The three least-certain parameters — `k_flip`,
   `β_leak`, `K_I` — can shift the leak/sensitivity balance; quantifying that shift is exactly
   what Module 7 exists to do. The conclusion should be stated as "at literature-anchored
   baselines," not as unconditional.

---

*All source values above were extracted from the fetched full-text PDFs, not from memory.*
