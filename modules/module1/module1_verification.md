# Module 1 — Verification Memo

**BioChronos** · Module 1 (Copper Sensor Dose-Response) · verification pass

## 1. Citation verification (CrossRef, this session)

All DOIs resolved against CrossRef; author lists, volume, pages, and article
numbers taken from CrossRef metadata. Corrections applied to the spec:

| Citation | Field | Draft (wrong) | CrossRef (correct) |
|---|---|---|---|
| Pang 2020 (WMC-007) | 1st author | "Chen" | **Pang** (Yilin Pang) |
| Pang 2020 | co-authors | "Wang, Fu" | **Ren, Li** (Pang, Ren, Li et al.) |
| Zahid 2012 (cus operon) | 1st author | "Zaheer" | **Zahid** (Nageena Zahid) |
| Zahid 2012 | co-authors | "Amin" | **Zulfiqar & Shakoori** |
| Fu 2024 (TCS biosensor) | article no. | "1373" (fabricated) | **1407** |

Verified as originally written: Stoyanov, Hobman & Brown 2001, *Mol. Microbiol.*
**39**:502–512; Fu et al. 2024, *Commun. Biol.* **7** (2024).

## 2. Source-access provenance (what was actually read)

| Source | Access | Read |
|---|---|---|
| Pang 2020 (Front. Microbiol. 10:3031) | Open (gold) | **Full text** — LOD 0.25 µM, linear range 0.39–78.68 µM extracted from text |
| Fu 2024 (Commun. Biol. 7:1407) | Open (gold) | **Full text** — WHO 30 µM / China 15 µM anchors, CusRS sensitivity |
| Stoyanov 2001 (Mol. Microbiol. 39:502) | Closed | **Abstract only** — proportionality claim taken from abstract |
| Zahid 2012 (Gene 495:81) | Closed | **Metadata only** — no full text/abstract retrieved |
| iScience 2025 (28:113715) | Landing only | **Abstract only** — cue/cus context |

## 3. Claims audit — grounded vs. withdrawn

- **GROUNDED (full text):** LOD 0.25 µM and linear range 0.39–78.68 µM (Pang 2020);
  regulatory anchors WHO 30 µM / China 15 µM (Fu 2024).
- **GROUNDED (abstract):** copA::lacZ "approximately proportional to the
  concentration of cupric ions in the medium" (Stoyanov 2001 abstract, verbatim) —
  supports non-cooperative n≈1.
- **GROUNDED (mechanism):** MerR-family regulators activate from a single
  metal-binding dyad ⟹ non-cooperative. (Textbook mechanism; Stoyanov confirms
  CueR is MerR-family.)
- **WITHDRAWN:** an earlier draft attributed a specific "apparent Hill coefficient
  1.0" to Zahid 2012. That paper's full text/abstract could not be retrieved
  (closed access), so the numerical figure is removed. The n≈1 baseline now rests
  on the Stoyanov abstract + MerR mechanism, and is swept 1–2 in Module 7 anyway.

## 4. Model self-consistency checks

- **Cross-validation:** combined-design 1 h net detection threshold = 0.264 µM,
  independently reproduces the measured WMC-007 LOD of 0.25 µM (Pang 2020) —
  the recorder model was built from Hsiao flip kinetics with no fitting to copper
  data, so this agreement is a genuine consistency check, not a circular fit.
- **Leak coupling:** leak-only (Cu=0) flip fractions match Module 2 exactly
  (untagged 11%/58%/99.9% at 1/4/24 h; combined 1%/9%/48%), confirming Models
  1 and 2 share one engine.
- **Model file** runs standalone and reproduces the results table.

## 5. Known discrepancy carried forward

Models 2 and 3 were built with a placeholder **n = 2**; Module 1 grounds the copper
sensor at **n ≈ 1**. Module 2's conclusions are n-independent (no-signal analysis),
but Model 3's dose-response / minimum-pulse panels used n = 2 and are technically
superseded. Logged as task T2 in OPEN_TASKS.md; Module 7's n-sweep (1–2) spans
both values.
