# BioChronos Dry Lab: Model-to-Wet-Lab Mapping Table

**As-of 2026-08-02.** Supersedes the device-to-model mapping in Parts 2-3 of the consolidated report, per M2 v3b Revision Instructions (Sections 2 and 3).

## Purpose
The consolidated report's model labels ("untagged", "weak RBS", "combined") do not cleanly identify the constructs the wet lab is now building. This table is the single authoritative crosswalk from each old model label to the current construct or experiment it maps to, and states explicitly what is a validated measurement versus a model-conditional prediction.

## 1. Naming convention (Section 3)
Computational sections are now **Dry Model D1, D2, ...**, not "Module 1-7". Wet-lab work uses the fixed names below.

| name | scope | construct | integrase / inversion |
|---|---|---|---|
| Wet M1 | copper sensor characterization | PCu-B0032-sfGFP | none |
| Wet M2 | pBAD memory test | pBAD/AraC-B0033-Bxb1-inversion reporter | Bxb1, arabinose-induced |
| Wet M3 | future integrated copper recorder | copper promoter driving Bxb1 | not yet finalized |
| Wet EXP-1 | (preserved wet-lab experiment name) | - | - |
| Wet EXP-2 | side-by-side M2 v3b-Untagged vs v3b-LVA | the two constructs below | Bxb1 |

## 2. The two current constructs (Section 1)
| construct | RBS | degradation tag | size | role |
|---|---|---|---|---|
| M2 v3b-Untagged | B0033 weak RBS | none | 5,871 bp | OFF-state comparison arm |
| M2 v3b-LVA | B0033 weak RBS | in-frame Bxb1 C-terminal AANDENYALVA | 5,904 bp | OFF-state comparison arm |

The primary decision comparison is **weak-RBS untagged vs the same weak RBS plus LVA**, differing only in the Bxb1 degradation term. Expected ON records are in-silico references, not order candidates.

## 3. Old model label -> current construct crosswalk
| consolidated-report label | maps to | status of the mapping | notes |
|---|---|---|---|
| "combined" (weak RBS + LVA) | **provisional proxy for M2 v3b-LVA** | proxy only | same RBS + LVA degradation term; leak parameter must be re-grounded to pBAD, not copper (see row below) |
| "weak RBS" (weak RBS, no tag) | **provisional proxy for M2 v3b-Untagged** | proxy only | same RBS, no tag |
| "untagged" (strong RBS, no tag) | **strong-RBS untagged model** | NOT a current construct | rename only; it is not M2 v3b-Untagged (that is weak-RBS). Keep as a modelling reference point, do not attribute to a built device |
| "LVA tag" (strong RBS + LVA) | strong-RBS + LVA model | NOT a current construct | modelling reference point only |

The re-mapped, pBAD-grounded comparison of the two current proxies is Dry Model D2 (delivered separately: D2_remap_comparison.md, D2_remap_predictions.png, D2_remap_model.py, D2_remap_results.csv).

## 4. Leak-parameter grounding by experiment (Section 3)
| experiment | correct basal-leak grounding | do NOT use |
|---|---|---|
| Wet M2 (current) | pBAD/AraC-specific leak (Guzman 1995 induced/repressed ratio; conditions = glucose repression / no inducer / arabinose) | copper-promoter leak (Fu 2024, beta_leak~100) |
| Wet M3 (future) | copper-promoter leak becomes relevant here | pBAD leak |

For Wet M2, "uninduced" means **no arabinose**, not no copper. Adding copper to M2 is not an informative induction test.

## 5. Status of the headline numbers from the old report
Every performance number below is a **model-conditional prediction**, not measured performance of the current sequences. None is validated for M2 v3b-Untagged or v3b-LVA until Wet EXP-2 data exist.

| headline number | old-report source | current status |
|---|---|---|
| 0.25 uM detection at 2 h | Module 1 (copper sensor) | prediction; belongs to the future copper recorder (Wet M3 / Wet M1 sensor), not M2 |
| 15.6 min to 1% false positive | Module 2 (combined) | prediction for the v3b-LVA proxy under a copper-leak assumption; must be recomputed with pBAD leak for M2 |
| 0.67 h trust horizon | Module 3 (combined) | prediction; copper-recorder context, not M2 |
| 94% three-state accuracy | Module 5 (combined single integrase) | prediction; multi-state copper readout, not M2 |

## 6. What Wet EXP-2 will actually measure (Section 5)
| output | reported as |
|---|---|
| GFP/OD600 and OD600 growth curves | time series per variant per condition |
| Orientation-specific endpoint colony PCR | ON-sequence detection frequency + state categories (OFF only / mixed / ON only) |

Endpoint colony PCR does NOT yield a molecular OFF/ON ratio or per-cell flip rate. A quantitative molecular fraction requires orientation-specific qPCR or ddPCR. After time-series data exist, both variants are fit simultaneously with experimental points, prediction intervals, and residuals, and the analysis states whether the data can actually distinguish the variants.
