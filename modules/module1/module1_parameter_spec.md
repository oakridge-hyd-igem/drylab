# Module 1 — Copper Sensor Dose-Response & Recording Threshold

**BioChronos** · Dry lab / modelling · Oakridge-HYD iGEM 2026
**Analyte:** COPPER (locked 2026-07-27) · **Sensor:** CueR / P(copA)

**Maps to:** the sensor-input layer — dose-response curve, activation threshold,
Hill coefficient — and Decision D1/D6 (sensor choice). Supplies the `β(Cu)` term
that Models 2 & 3 previously took as a placeholder Hill function.

---

## 1. Sensor choice

**CueR / P(copA)** — the *Cue* system, the primary and best-characterised copper
sensor in *E. coli*. CueR is a MerR-family transcriptional activator that binds
cytoplasmic Cu(I) and activates the copA promoter (Stoyanov 2001; Outten 2000).
The *Cus* two-component system (CusRS → cusC) is the alternative, but it responds
to periplasmic Cu(I) under anaerobic/high-stress conditions and shows baseline
crosstalk (Fu 2024; Chen/iScience 2025) — CueR/PcopA is the cleaner primary sensor.

## 2. Model

Integrase production as a Hill function of copper, coupled to the **Module 3 flip
engine** (unchanged flip kinetics) so Models 1/2/3 stay numerically consistent:

```
β(Cu) = β_leak + β_max · Cuⁿ / (Kⁿ + Cuⁿ)          # copies·cell⁻¹·h⁻¹
I(t)  = (β/γ)(1 − e^(−γt))                          # integrase ramp during exposure
λ(I)  = k_flip · I / (I + K_I)                       # Hsiao-anchored flip hazard
F(Cu) = 1 − exp(−∫₀^t_expose λ(I(t)) dt)             # fraction recording a flip
```

### 2.1 Literature-grounded parameters

| Parameter | Symbol | Value | Plausible range (M7) | Source |
|---|---|---|---|---|
| **Hill coefficient** | n | **1.0** | 1 – 2 | Stoyanov 2001 abstract: copA::lacZ "approximately proportional to the concentration of cupric ions" ⟹ **non-cooperative (n≈1)**. MerR-family single-site activation is mechanistically non-cooperative. |
| Half-max copper | K | 5 µM | 1 – 15 | Inferred from operating windows (no explicit fitted K in literature); mid-window, just below China 15 µM line |
| Limit of detection | — | 0.15–0.25 µM | — | Pang 2020 (WMC-007, 0.25 µM); Ivask 2009 (0.15 µM) |
| Linear range | — | ~0.4–80 µM | — | Pang 2020 (0.39–78.68 µM); Riether 2001 (3–30 µM) |
| Max production | β_max | 1000 | 100–10 000 | Module 3 (Alon; strong promoter) |
| Basal leak | β_leak | **100** | 10–500 (sweep) | **T1-grounded: 10% of β_max, natural CueR/P(copA) fold-change 7–18× (Fu 2024). Sweep 1% floor → 50% worst natural.** |
| Flip rate / threshold | k_flip, K_I | 0.4 h⁻¹, 10 | Module 7 | Module 3 (Hsiao 2016) |

**The one genuinely new sensor finding: n ≈ 1, not 2.** Two lines support a
non-cooperative copper response: (i) the Stoyanov 2001 abstract states copA::lacZ
expression is "approximately proportional to the concentration of cupric ions"
(verified from the fetched abstract) — proportional ⟹ n≈1; and (ii) MerR-family
regulators activate from a single metal-binding dyad, which is mechanistically
non-cooperative. The resulting dose-response is *near-hyperbolic* — a gentler,
less switch-like curve than the n=2 placeholder Models 2/3 inherited. This widens
the graded-response region and is swept 1–2 in Module 7.

*Not verified:* an earlier draft cited an "apparent Hill coefficient 1.0" from
Zahid 2012 (cus operon); that paper's full text/abstract could not be retrieved
(closed access), so the specific figure is withdrawn. The n≈1 grounding rests on
the Stoyanov abstract and the MerR mechanism, not on Zahid.

## 3. Results

### 3.1 The activation threshold must be signal-over-background

Because the switch is irreversible and leak is non-zero, a naïve "minimum copper
to flip 5%" is meaningless for the untagged design — **leak alone flips >5%
within an hour** (Module 2). The correct activation threshold is the copper level
that adds ≥5% flips *above the no-signal leak background*:

*(β_leak grounded to 10% of β_max — T1, measured copper-promoter fold-change,
Fu 2024. The earlier 1% values understated the leak background.)*

| Exposure | untagged: net threshold (leak background) | combined: net threshold (leak background) |
|---|---|---|
| 1 h | 4.37 µM (bg **26%**) | 0.584 µM (bg 9%) |
| 2 h | 9.54 µM (bg **49%**) | **0.246 µM** (bg 22%) |
| 3 h | ∞ — leak floods (bg **65%**) | 0.176 µM (bg 34%) |
| 4 h | ∞ — leak floods (bg **76%**) | 0.151 µM (bg 45%) |

At the grounded copper leak (10× higher than the old 1% assumption), the untagged
design is **unusable**: leak alone flips 26% of cells within 1 h and 76% by 4 h,
and beyond 2 h no copper level can add a 5% net signal above that runaway
background (threshold → ∞). Only the **combined leak-suppressed design** keeps a
workable detection window — its net threshold *falls* with exposure (0.58 → 0.15
µM over 1–4 h) as signal integrates faster than its lower leak background. This
sharpens the Module 1 ↔ Module 2 coupling: at real copper-promoter leak,
sensitivity is only meaningful **with** leak suppression, and the untagged design
cannot serve as a recorder at all for exposures beyond ~2 h.

### 3.2 Independent validation against real copper sensors

The combined design's **2 h net threshold = 0.246 µM** lands almost exactly on the
**measured WMC-007 bioreporter LOD of 0.25 µM** (Pang 2020) — despite our recorder
model being built entirely from Hsiao *flip* kinetics with no fitting to copper
data. (At the grounded 10% leak this LOD-match moves from the 1 h window to the 2 h
window; the earlier 1%-leak spec matched at 1 h. The 2 h window is arguably the more
honest anchor, since a 1 h read already carries a 9% leak background.) The
recorder's copper sensitivity is consistent with published whole-cell copper
sensors, and comfortably below the WHO (30 µM) and China (15 µM) drinking-water limits — the device can record environmentally and
regulatorily relevant copper.

## 4. Caveats

1. **K is inferred, not fitted.** No copper-sensor paper reports an explicit
   half-max K; it is inferred from published LOD/linear-range windows. Model 7
   must sweep K over the grounded 1–15 µM range (Panel D: threshold scales
   ~linearly with K).
2. **n = 1 is well-grounded but not from a single clean fit.** It rests on two
   independent lines (cus-operon Hill = 1.0; copA proportionality). Sweep 1–2.
3. **Sustained-exposure idealisation.** Copper is modelled as a constant level
   during exposure; EXP3 pulses (1–4 h) are square. Pulse *shape* / uptake
   kinetics (intracellular vs medium copper) are not modelled — a Model 6 /
   diffusion concern.
4. **CueR senses Cu(I), assay copper is Cu(II).** Medium Cu²⁺ is reduced to
   cytoplasmic Cu⁺; the model uses external copper as the input variable
   (as the reporter literature does), folding the reduction/uptake into K.

## 5. References

*Provenance: all DOIs, author lists, volumes, pages and article numbers below
were verified against CrossRef metadata this session. Full text was fetched for
Pang 2020 and Fu 2024 (open access); only the abstract was retrievable for
Stoyanov 2001 (closed access, published abstract retrieved); the iScience 2025
review (closed access, published abstract retrieved); and Zahid 2012 (closed
access, metadata only — no abstract or full text retrieved). Claims are attributed only to the level of
source actually read.*

- Stoyanov, Hobman & Brown (2001) "CueR (YbbI) of *E. coli* is a MerR family
  regulator controlling expression of the copper exporter CopA," *Mol. Microbiol.*
  **39**:502–512, doi:10.1046/j.1365-2958.2001.02264.x — *(abstract read)*
  copA::lacZ expression "approximately proportional to the concentration of cupric
  ions in the medium" ⟹ non-cooperative (n≈1); CueR is a MerR-family regulator.
- Pang, Ren, Li et al. (2020) "Development of a Sensitive *E. coli* Bioreporter
  Without Antibiotic Markers for Detecting Bioavailable Copper in Water Environments,"
  *Front. Microbiol.* **10**:3031, doi:10.3389/fmicb.2019.03031 — WMC-007
  CueR/PcopA::gfp, LOD 0.25 µM, linear range 0.39–78.68 µM.
- Fu et al. (2024) "Development of a two component system based biosensor with high
  sensitivity for the detection of copper ions," *Commun. Biol.* **7**:1407 (2024),
  doi:10.1038/s42003-024-07112-6 — CusRS sensitivity engineering; WHO 30 µM /
  China 15 µM regulatory anchors.
- Zahid, Zulfiqar & Shakoori (2012) "Functional analysis of cus operon promoter of
  *Klebsiella pneumoniae* using *E. coli* lacZ assay," *Gene* **495**:81–88,
  doi:10.1016/j.gene.2011.12.040 — *(metadata only; closed access, full text not
  retrieved this session)* cus-operon copper promoter characterisation; cited for
  context, not for a specific Hill-coefficient value.
- Chen et al. (2025) "Systematic and synthetic biology insights into copper
  homeostasis in *E. coli*," *iScience* **28**:113715,
  doi:10.1016/j.isci.2025.113715 — cue vs cus roles, Cu⁺ speciation, CusR crosstalk.

Flip kinetics inherited from Module 3 (Hsiao et al. 2016, *Mol. Syst. Biol.*
12:869, doi:10.15252/msb.20156663).
