# BioChronos Dry Lab: Parameter Provenance Table

**As-of 2026-08-02.** Per M2 v3b Revision Instructions Section 5. Every model parameter is labelled as one of **measured / literature proxy / assumed / fitted**, with source, units, host strain / construct context, medium and temperature, and applicable time scale. Unknown current values remain TBD.

## Status of the four classes right now
- **measured:** none yet. No wet-lab measurement of any BioChronos construct exists. These fill in after Wet EXP-1 (sensor), Wet EXP-2 (pBAD memory), and later Wet M3.
- **fitted:** none yet. Fitting requires Wet EXP-2 time-series data; both M2 variants will be fit simultaneously once data exist.
- **literature proxy:** values taken from published work on related (not identical) systems, used as stand-ins until measured.
- **assumed:** engineering estimates or scale-setting constants with no direct literature value; these carry the most model risk and are the first targets for wet-lab replacement.

All current parameters are therefore **literature proxy or assumed** (13 proxy, 5 assumed). No headline prediction should be read as a measured property of a BioChronos construct.

## Parameter table

| parameter | value | units | class | source | strain / context | medium, temp | time scale | notes |
|---|---|---|---|---|---|---|---|---|
| k_flip (Bxb1 flip rate) | 0.4 | 1/h | literature proxy | Hsiao et al. 2016 Mol Syst Biol 12:869 (nominal fit) | E. coli, Bxb1 integrase-memory construct | not stated / 37 C typical | per-hour hazard | half-saturating hazard k_flip*I/(I+K_I) |
| K_I (integrase half-saturation) | 10 | molecules/cell | literature proxy | Hsiao et al. 2016 (Kd ~10 molecules) | E. coli | 37 C | steady-state | first-order proxy for tetramer assembly |
| gamma_dil (dilution) | 0.30 | 1/h | literature proxy | doubling time ~2.3 h -> ln2/2.3 | E. coli exponential growth | rich, 37 C | per generation | growth-rate dilution of integrase |
| gamma_LVA (ssrA-LVA degradation) | 1.04 | 1/h | literature proxy | Andersen et al. 1998 AEM 64:2240 (ssrA-LVA t1/2 ~40 min) | E. coli | 37 C | per-hour | adds to dilution for tagged variant (total 1.34) |
| beta_max (strong-RBS induced integrase) | 1000 | copies/cell/h | assumed (scale-setting) | normalization constant | E. coli | - | - | reference production scale; only ratios matter |
| f_rbs (B0033 weak RBS factor) | 0.1 | dimensionless | assumed | B0033 ~10x weaker than strong RBS (order-of-magnitude) | B0033 weak RBS | - | - | scales induced AND basal equally |
| beta_leak copper (P(copA)/CueR) | ~100 (10% of beta_max) | copies/cell/h | literature proxy | Fu et al. 2024 Commun Biol 7:1407 (natural CueR fold-change 7-18x) | CueR/P(copA), E. coli | not stated | steady-state | FUTURE Wet M3 only; NOT for current M2 |
| pBAD leak (glucose repression) | 1/1200 (0.083%) | fraction of induced | literature proxy | Guzman et al. 1995 J Bacteriol 177:4121 (induced/repressed up to 1200x) | pBAD/AraC, E. coli | glucose, 37 C | steady-state | current Wet M2 OFF-state (tightest) |
| pBAD leak (no inducer) | ~0.5% (0.1-2%) | fraction of induced | assumed | bracket; no arabinose, no glucose | pBAD/AraC | 37 C | steady-state | current Wet M2; EXP-2 uninduced fits this |
| K_copper (sensor half-max) | ~5 (1-15) | uM Cu(II) | assumed | inferred from operating windows; no fitted literature value | CueR/P(copA) | not stated | dose-response | sensor Hill K; Wet M1 dose-response fits it |
| n_hill (sensor Hill coefficient) | 1.0 | dimensionless | literature proxy | Stoyanov et al. 2001 Mol Microbiol 39:502 (copA proportional to Cu; MerR single-site) | CueR/P(copA) | not stated | dose-response | MerR-family single-site mechanism |
| copper LOD (sensor) | 0.25 | uM Cu(II) | literature proxy | Pang et al. 2020 Front Microbiol 10:3031 (WMC-007) | chromosomal copAp::gfpmut2, delta-copA delta-cueO delta-cusA | LB, 5 h, 37 C | 5 h assay | sensor detection limit, NOT DH5a host tolerance |
| fitness cost c | 2% (0-5%) | growth-rate fraction | literature proxy | Bienick 2014 PLoS ONE 9:e109105 (linear expression-growth); Gonzalez-Colell 2025 Sci Rep 15:31311 (Bxb1 burden ~0) | E. coli | glucose minimal / rich | per generation | Module 4 memory-stability selection term |
| extrinsic noise CV | 0.3 | dimensionless | literature proxy | Taniguchi et al. 2010 Science 329:533 (extrinsic floor eta_p^2 >= 0.1) | E. coli proteome | not stated | single-cell | Module 5 multi-state separability; T3 wants EXP3 confirmation |
| F_detect (readout threshold) | 0.05 | flipped fraction | assumed | 5% population flip = detectable | - | - | readout | trust-horizon / detection threshold definition |
| Cu host-tolerance onset (WT) | >8 | uM Cu(II) | literature proxy | Macomber & Imlay 2009 PNAS 106:8344 (Fig 1A) | WT W3110 | simple salts + glucose, 37 C | growth assay | host growth-tolerance, NOT sensor range |
| Cu host MIC (WT, rich) | 3.0 | mM CuSO4 | literature proxy | Franke et al. 2003 J Bacteriol 185:3804 (GR1 delta-cueO control) | W3110 delta-cueO | optimized LB pH7.5, 37 C | MIC | conservative lower bound; DH5a cueO+ tolerates >=3 mM |
| Cu host MIC (WT, minimal) | 8 | mM CuSO4 (48 h) | literature proxy | Rosenberg et al. 2025/2026 (biorxiv 2025.08.27.672559) | WT BW25113 | MOPS minimal, 48 h, 37 C | MIC | minimal-medium MIC; chelation-dependent |

## Isolated-variable rule for the LVA comparison (Section 5)
For the M2 v3b-Untagged vs v3b-LVA comparison, **the only parameter allowed to differ is the effective Bxb1 degradation rate** (gamma: 0.30/h untagged vs 1.34/h LVA). Promoter strength, RBS strength (both B0033, f_rbs=0.1), copy number, and flipping kinetics (k_flip, K_I) are held identical. This is confirmed in code (D2_remap_model.py).

## Wet EXP-2 outputs and how they update this table
| output | fills in |
|---|---|
| GFP/OD600 and OD600 growth curves | fitness cost c (measured); induced expression scale |
| Orientation-specific endpoint colony PCR (ON-frequency + state categories: OFF only / mixed / ON only) | qualitative OFF/ON separation per variant |
| (if quantitative fraction needed) orientation-specific qPCR or ddPCR | molecular flip fraction -> fits k_flip, pBAD leak fractions |

Endpoint colony PCR does not yield a molecular OFF/ON ratio or per-cell flip rate. After time-series data are available, both variants are fit simultaneously (points, prediction intervals, residuals), and the analysis reports whether the data can actually distinguish the variants. The `class` column for the fitted parameters moves from "literature proxy"/"assumed" to "fitted" at that point.

## Highest-risk assumed parameters (replace first)
1. **pBAD leak (no inducer) ~0.5%** - assumed bracket; directly sets the OFF-state false-positive prediction for Wet M2. Wet EXP-2 uninduced data fits it.
2. **K_copper ~5 uM** - assumed sensor half-max; no fitted literature value. Wet M1 dose-response fits it.
3. **f_rbs 0.1 and beta_max 1000** - order-of-magnitude scale constants; only their ratio enters, but the ratio sets the induced-vs-basal scale.
