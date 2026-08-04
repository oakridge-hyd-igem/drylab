# BioChronos Dry Lab: Deliverables Index for Prof. Nishida and the Wet Lab

**As-of 2026-08-02.** This index maps each item the Dry Lab was asked for (M2 v3b Revision Instructions, Section 10) to the exact artifact that satisfies it, with a status. The consolidated report is the entry point; the rest are its companion deliverables.

## Status of the Section 10 required deliverables

| # | Requested deliverable (Notion Section 10) | Status | Artifact |
|---|---|---|---|
| 1 | One-page comparison of M2 v3b-Untagged and M2 v3b-LVA | Delivered | D2_remap_comparison.md (+ model .py, results .csv) |
| 2 | Corrected model-to-wet-lab mapping table | Delivered | D2_model_wetlab_mapping.md |
| 3 | Parameter table with provenance, units, construct context, status | Delivered (all params literature-proxy/assumed; measured/fitted fill in after wet data) | parameter_provenance.md / .csv |
| 4 | Prediction plots, 3 induction conditions x 2 variants | Delivered | D2_remap_predictions.png |
| 5 | Numerical-audit appendix (Section 6 issues A, B, C) | Delivered | section6_audit_appendix.md (+ figure) |
| 6 | Citation-verification log for the F4 prior-art claims | Delivered | F4_citation_verification_log.md |
| 7 | Conclusion table (supported / inconclusive / contradicted) | Delivered | claim_status_rollup.md |
| 8 | One next-iteration proposal, only after Wet EXP-2 data | Not yet due (gated on Wet EXP-2 data, per the instruction itself) | pending Wet EXP-2 |

Seven of eight are delivered. Item 8 is intentionally open: the instruction gates it on Wet EXP-2 data, which do not exist yet.

## Also delivered (Section 4, DNA-free wet-lab week)
| Deliverable | Status | Artifact |
|---|---|---|
| CuSO4 growth-tolerance ladder + plate-reader settings + plate layout + model linkage | Delivered | dnafree_week_deliverables.md |
| 96-well plate layout figure | Delivered | cu_plate_layout.png |
| Grounded copper-tolerance anchors (data file) | Delivered | cu_tolerance_grounded.json |

## Contents of this package

**Primary**
- consolidated report (Revision 2) - the entry point
- D2_model_wetlab_mapping.md - how old model labels map to the current constructs
- claim_status_rollup.md - what is supported vs pending vs corrected

**Dry Model D2 (the current M2 comparison)**
- D2_remap_comparison.md, D2_remap_predictions.png, D2_remap_model.py, D2_remap_results.csv

**Provenance and audit**
- parameter_provenance.md / .csv, section6_audit_appendix.md, section6_audit.png, F4_citation_verification_log.md

**DNA-free wet-lab week**
- dnafree_week_deliverables.md, cu_plate_layout.png, cu_tolerance_grounded.json

## Read-before-use notes
- Every headline number (0.25 uM, 15.6 min, 0.67 h, 94%) is a model-conditional prediction, not measured performance of the current sequences.
- The current M2 comparison is a pBAD/AraC memory test; its leak parameter is pBAD-specific, not the copper-promoter value.
- Section 10 item 8 (next-iteration proposal) is delivered only after Wet EXP-2 data exist.
- Not part of Section 10, noted for completeness: Module 6 (deployment/encapsulation) is not built; the Gold Medal criterion mapping (Q1) and novelty determination (Q2) are team/registration-owned.
