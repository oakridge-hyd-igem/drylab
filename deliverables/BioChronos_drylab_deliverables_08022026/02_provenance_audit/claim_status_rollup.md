# BioChronos Dry Lab: Project-Wide Claim Status Roll-Up

**As-of 2026-08-02.** Per M2 v3b Revision Instructions Section 7. Every substantive dry-lab claim is classified as:
- **SUPPORTED** - confirmed by published data, an internal numerical audit, or a code-level check this project performed.
- **INCONCLUSIVE / PENDING** - a model-conditional prediction with no wet-lab test yet; direction plausible but unvalidated.
- **CONTRADICTED / CORRECTED** - found wrong (by audit, re-grounding, or citation check) and since corrected.

No claim about a BioChronos construct's performance is yet SUPPORTED by a wet-lab measurement; "SUPPORTED" here means supported by literature/audit/code, not by our own bench data.

## SUPPORTED (literature / audit / code)
| # | claim | basis |
|---|---|---|
| S1 | LVA tag reduces steady-state integrase level and OFF-state leak vs untagged | degradation adds to dilution (gamma 0.30 -> 1.34/h); Andersen 1998 ssrA-LVA half-life; confirmed in D2 code |
| S2 | For the M2 v3b comparison, the two variants differ only in the degradation term | code-verified (D2_remap_model.py): identical f_rbs, k_flip, K_I, gamma_dil; LVA adds gamma_LVA only |
| S3 | Time-to-threshold scales ~1/sqrt(beta_leak), not linearly with steady-state leak | Section 6 numerical audit A (coarse-vs-fine reproduced within rounding; t proportional to 1/sqrt(beta_leak)) |
| S4 | Leak is the single dominant lever on the trust window and detection threshold | Module 7 OAT sensitivity: beta_leak swings A_trust 231%, B_thresh 1540%; all other params far smaller |
| S5 | Copper sensor LOD ~0.25 uM is achievable and below regulatory limits (WHO 30, China 15 uM) | Pang 2019 WMC-007 LOD 0.25 uM (verified full text); regulatory limits from standards |
| S6 | Host copper tolerance is medium-dependent and set by chelation, spanning ~1000-fold | Franke 2003 (3 mM LB), Macomber 2009 (>8 uM simple-salts, verified PDF), Rosenberg (8 mM MOPS) |
| S7 | Base part BBa_K907000 is KAIST 2012 (Trc/IPTG), not Peking 2012 | Registry XML (entered 2012-09-20, Choe & Lee); F4 citation log |
| S8 | Memory (inverted DNA) is heritable and stable over ~90 generations under modelled cost | Module 4: 90-gen retention ~100% with leak on, ~60% leak-free worst case (audit C reconciled) |

## INCONCLUSIVE / PENDING (model-conditional, needs wet lab)
| # | claim | test that would resolve it |
|---|---|---|
| P1 | LVA improves net OFF/ON separation (cuts OFF leak ~3x at ~0.2 pt ON cost) | Wet EXP-2 side-by-side GFP/OD + colony PCR |
| P2 | 15.6 min to 1% false positive (v3b-LVA proxy) | Wet EXP-2, AND recompute with pBAD (not copper) leak first |
| P3 | 0.67 h trust horizon; sub-minute pulses recordable | future copper-recorder wet test (Wet M3) |
| P4 | 94% three-state accuracy at 1.25 h readout; 2nd integrase adds ~4 pts | multi-state copper readout experiment (not M2); needs EXP3 noise CV |
| P5 | pBAD OFF-state leak ~0.5% (no inducer), ~0.08% (glucose) | Wet EXP-2 uninduced arm |
| P6 | K_copper ~5 uM sensor half-max | Wet M1 copper dose-response |
| P7 | Fudan China 2017 owns the pBAD/B0033/AraC integrase characterization | manual Registry page + team-wiki confirmation (HTML page was 403 this session) |
| P8 | Waseda 2020 ssrA-LVA; BBa_K2384014 multi-repressor; TP901 lowest-leak | manual Registry/wiki checks (F4 log C2, C4, C5) |
| P9 | Single leak-suppressed integrase sufficient for a binary recorder | wet recorder test; part-level suppression alone is ~5000x short of week-long leak-tight (Module 2) |

## CONTRADICTED / CORRECTED (found wrong, fixed)
| # | original claim | correction | where |
|---|---|---|---|
| C1 | Current M2 leak grounded in copper-promoter value (beta_leak~100) | M2 is a pBAD/AraC test; leak must be pBAD-specific (Guzman 1995). Copper leak reserved for future Wet M3 | Section 3 / D2 re-map |
| C2 | Report's "combined"/"weak RBS"/"untagged" map directly to built devices | only proxies; "untagged" was strong-RBS, renamed strong-RBS untagged model, not a current construct | mapping table |
| C3 | "44x lower steady-state leak -> 44x longer time to threshold" | Iss is 44x lower but hazard lambda only ~2.3x lower; time scales ~1/sqrt(beta_leak) | Section 6 audit A |
| C4 | BBa_K907000 = Peking 2012, pBAD/B0033 | KAIST 2012, Trc/IPTG; pBAD/B0033 is later work (Fudan 2017 per instructions, team linkage TBD) | F4 log C1 |
| C5 | CFLS 2025 LAA-LAA supports Bxb1 degradation-tagging | CFLS 2025 context is NhaA/toxin switch, not integrase; general prior art only | F4 log C3 |
| C6 | Macomber anchors "cited, not verified" / MOPS medium / LEM33 a fabrication | verified from PDF: WT W3110 >8 uM simple-salts glucose (not MOPS); LEM33 is the real triple-KO strain | DNA-free deliverables |
| C7 | "No prior iGEM team quantified cumulative leak" (absolute) | not defensible without a documented systematic search; reframed to "in pages reviewed here" | F4 log C6 |
| C8 | module 5 combined beta_max = 200 | corrected to 100 (weak RBS x0.1 scales induced and basal equally) | Module 5 |

## One-line project status
The dry lab has a coherent, literature-grounded model whose central lesson (leak dominates; LVA helps; memory is stable) is SUPPORTED by literature/audit/code, but **every quantitative performance claim about a BioChronos construct is PENDING wet-lab data.** The main corrections since the first report were re-grounding M2 leak to pBAD (not copper), re-mapping model labels to the actual constructs, and fixing citation attributions.
