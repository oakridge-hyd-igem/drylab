# BioChronos Dry Model D2 (re-mapped): M2 v3b-Untagged vs M2 v3b-LVA

**As-of 2026-08-02. Section 2 re-mapping per M2 v3b Revision Instructions.**

## What changed from the old report

The old report compared **strong-RBS untagged** against **combined** (weak RBS + LVA), which confounds two design changes. The current wet-lab decision (Wet EXP-2) is a controlled comparison of two constructs differing in **one** feature:

- **M2 v3b-Untagged:** B0033 weak RBS, no degradation tag (5,871 bp)
- **M2 v3b-LVA:** identical plus an in-frame Bxb1 C-terminal AANDENYALVA tag (5,904 bp)

**Isolated-variable confirmation (in code).** The two model conditions differ only in the LVA degradation term: both use f_rbs=0.1, and share k_flip=0.4/h, K_I=10, gamma_dil=0.3/h. The only difference is gamma_LVA (+1.04/h) on the LVA variant. No promoter, RBS, copy-number, or kinetic parameter differs. This is the correct weak-RBS-untagged vs weak-RBS+LVA pairing.

## Leak is pBAD-grounded, NOT copper

Current M2 is a **pBAD/AraC arabinose** memory test, so the basal-leak parameter is pBAD-specific. The copper-promoter leak (Fu 2024, beta_leak~100, ~10% of induced) is **reserved for the future integrated copper recorder (Wet M3)** and is not used here. pBAD leak fractions:

| condition | leak fraction | basis | provenance |
|---|---|---|---|
| glucose repression | 1/1200 (0.083%) | Guzman et al. 1995, J Bacteriol 177:4121-4130 (induced/repressed up to 1200-fold) | literature proxy |
| no inducer (no ara, no glucose) | ~0.5% (bracket 0.1-2%) | leakier OFF without catabolite repression | assumed (TBD, EXP-2 fits) |
| arabinose induction (ON) | 100% (defines beta_max) | full induction | defines scale |

Note "uninduced" for M2 means **no arabinose**, not no copper. Adding copper to M2 is not an informative induction test.

## Predictions (both variants x 3 conditions)

Flipped fraction (%), fine-grid model:

| condition | variant | 8 h | 16 h |
|---|---|---|---|
| OFF glucose repression | v3b-Untagged | 5.3 | 12.9 |
| OFF glucose repression | v3b-LVA | 1.8 | 3.7 |
| OFF no inducer | v3b-Untagged | 25.6 | 52.3 |
| OFF no inducer | v3b-LVA | 9.9 | 19.7 |
| ON arabinose | v3b-Untagged | 94.9 | 99.8 |
| ON arabinose | v3b-LVA | 93.6 | 99.6 |

## The design question the instructions raised: does LVA cost induced switching?

The revision instructions correctly warned: "Do not assume LVA is superior. It may reduce basal Bxb1 accumulation while also reducing induced switching." The model quantifies both sides:

- **OFF-state benefit (large):** LVA cuts leak-driven false positives roughly 3-fold. Glucose repression at 16 h: 12.9% -> 3.7%. No-inducer at 16 h: 52.3% -> 19.7%.
- **ON-state cost (small):** LVA lowers induced flip only marginally, because at full induction the integrase is far above the K_I=10 half-saturation, so the hazard is near-saturated even with the higher degradation. Arabinose ON at 16 h: 99.8% -> 99.6%; at 8 h: 94.9% -> 93.6%.

**Net (model prediction, to be tested):** for pBAD-driven Bxb1, the LVA tag improves OFF/ON separation because the induced integrase level sits on the saturated part of the flip-hazard curve, so trimming its half-life barely touches the ON signal while sharply cutting the OFF leak. This is a **model-conditional prediction**, not measured performance, and it hinges on the assumed no-inducer leak fraction and on whether the LVA tag's real effect on induced switching matches the gamma-only model. Wet EXP-2 time-series data are required to confirm it.

## Caveats
- All flip percentages are model-conditional predictions for these sequences, not measured values.
- The no-inducer leak fraction (0.5%) is assumed; EXP-2 uninduced data should fit it.
- The gamma-only model assumes LVA changes only the degradation rate. If the tag also perturbs folding, activity, or expression, the ON-state cost could be larger than predicted; EXP-2 tests this directly.
- Endpoint colony PCR gives state categories (OFF only / mixed / ON only), not a molecular flip rate; a quantitative fraction needs orientation-specific qPCR or ddPCR.
