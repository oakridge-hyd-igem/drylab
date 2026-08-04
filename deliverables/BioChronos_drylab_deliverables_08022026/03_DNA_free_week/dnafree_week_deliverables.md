# BioChronos Dry Lab: DNA-Free Wet-Lab Week Deliverables

**As-of: 2026-08-02. Prepared for the DNA-free wet-lab week beginning 2026-08-03.**
**Scope: untransformed DH5alpha copper growth-tolerance screen. No DNA, no integrase, no inversion.**

This package answers Section 4 of the M2 v3b Revision Instructions. It is a host growth-tolerance screen: it constrains the safe copper range for later Wet EXP-1 and Wet M3. It is NOT a measurement of Bxb1 flipping or pBAD leakage.

---

## 1. CuSO4 growth-tolerance ladder for untransformed DH5alpha

### Why not the recorder ladder
The 0.1-20 uM ladder proposed elsewhere is the **sensor operating window** (CueR/P(copA) dose-response), designed for inversion measurements. It is the wrong scale for host growth tolerance and must not be reused here. Host copper tolerance is set by the **medium**, and it sits well above the sensor range in rich media.

### Grounded anchors
| source | medium | strain | reported value | note |
|---|---|---|---|---|
| Franke, Grass, Rensing, Nies 2003, J Bacteriol 185:3804-3812 | optimized LB pH 7.5, aerobic | W3110 delta-cueO (control strain GR1) | MIC 3.0 mM CuSO4 | conservative lower bound (cueO- is more Cu-sensitive than WT; true DH5alpha cueO+ tolerates >=3 mM) |
| Rosenberg, Umerov, Teaer, Ivask 2025/2026, Microbiol Spectrum (biorxiv 2025.08.27.672559) | MOPS minimal medium (MMM), aerobic | WT BW25113 | MIC 8 mM Cu over 48 h; 4 mM = 0.5xMIC | mM-scale MIC in a MINIMAL medium (see chelation note). A "4 mM at 16-20 h" endpoint is sometimes quoted but is not in the source; do not rely on it |
| Macomber & Imlay 2009, PNAS 106:8344-8349 | simple salts medium, glucose sole carbon source, aerobic | WT W3110 | growth defect at Cu(II) >8 uM (Fig 1A; doses tested 0, 8, 16, 32 uM) | direct quote: "WT E. coli exhibited a growth defect when Cu(II) concentrations exceeded 8 uM" |
| Macomber & Imlay 2009, PNAS 106:8344-8349 | simple salts medium, glucose sole carbon source, aerobic | LEM33 (copA cueO cusCFBA triple KO) | growth defect at Cu(II) 0.25-1 uM (Fig 1B; doses tested 0, 0.25, 0.5, 1 uM) | context only, NOT our strain (copper-hypersensitive knockout) |

**Strain note.** DH5alpha carries wild-type copper homeostasis (copA, cueO, cus). It is bracketed by the WT anchors above, NOT by the triple-knockout 0.25 uM figure. The 0.25 uM value from the sensor literature is Pang et al.'s WMC-007 LOD, measured with a chromosomal copAp::gfpmut2 reporter in a delta-copA delta-cueO delta-cusA strain on a 5 h assay. That is a sensor detection limit in a copper-hypersensitive knockout background, not a host growth-tolerance figure for the WT-homeostasis DH5alpha, and it is used here only as external context.

Note on a coincidence: 0.25 uM appears twice above, once as the lowest Cu(II) dose at which Macomber's LEM33 triple-KO shows a growth defect and once as Pang's WMC-007 biosensor LOD. These are two independent values (a growth-assay dose vs a biosensor detection limit) that coincide because both use the same copper-hypersensitive delta-copA delta-cueO delta-cusA (cusCFBA) background; neither describes DH5alpha host tolerance.

### Proposed ladders (medium-dependent; run whichever medium the plate assay uses)
Two ladders are provided because copper host-tolerance is set by the **medium's copper-chelating / bioavailability capacity**, not simply by a rich-vs-minimal label. Rich media (LB) chelate copper heavily via amino acids and peptides, raising the tolerated dose into the mM range; a phosphate/MOPS glucose-minimal medium with no amino acids leaves copper far more bioavailable, so WT onset can be ~1000-fold lower (~8 uM, Macomber 2009). Caution: this is not a clean rich/minimal dichotomy. Rosenberg 2025/2026 measured an 8 mM MIC in a **MOPS minimal** medium, so the two minimal-medium anchors (Macomber's >8 uM onset in a simple-salts glucose medium vs Rosenberg's 8 mM MIC in MOPS minimal) span ~1000-fold depending on medium chemistry. Each ladder is designed to put the growth-inhibition onset **inside** the ladder, not at an edge.

- **LB (rich) ladder, mM CuSO4:** 0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0
  Fine steps bracket the 1-3 mM onset/MIC band (Franke conservative lower bound 3 mM MIC); the 8-10 mM top exceeds the 48 h MIC (8 mM, Rosenberg 2025/2026) so full inhibition is captured even at the longer read.
- **M9 (minimal) ladder, uM CuSO4:** 0, 1, 2, 5, 10, 25, 50, 100, 250, 500
  Log-spaced; brackets the WT onset (Macomber 2009: growth defect at Cu(II) >8 uM) up past the KO range. **De-risking (important):** the ~8 uM onset was measured in Macomber's "simple salts medium with glucose as sole carbon source," which is not necessarily identical to the wet lab's M9 formulation. Minimal-medium copper tolerance spans ~1000-fold across these anchors (Rosenberg's 8 mM MIC in MOPS minimal vs Macomber's >8 uM onset in a simple-salts glucose medium), because tolerance is set by the medium's copper bioavailability / chelation capacity, not the minimal label. If the wet lab's M9 chelates copper more than Macomber's medium did, the true onset could sit well above the ~8 uM anchor and the ladder could miss it. Before committing, run a coarse range-finder (0, 10, 100, 1000 uM) in the actual M9 formulation, or widen the ladder top. State the exact M9 composition on the bench sheet.

### Assumptions and conditions (state these on the bench sheet)
- Medium: LB or M9-glucose, as chosen by wet lab. Report pH (copper toxicity rises near pH 7.5).
- Exposure: 37 C, aerobic, shaking; endpoint OD600 at a fixed time (recommend both a 16-20 h read and a 48 h read). The reference MIC is Rosenberg's 8 mM over 48 h; reading at two timepoints captures the onset without assuming how the MIC shifts with incubation time.
- CuSO4 prepared fresh from a filter-sterilised stock; add to medium immediately before inoculation.
- **Safety / supervisor confirmation required** before handling millimolar copper stocks; confirm disposal route for copper-containing waste with Prof. Nishida / lab safety officer.

---

## 2. Plate-reader settings proposal

**Instrument model: TBD.** Confirm the exact plate reader before finalising gain and filter/monochromator settings. The values below are generic sfGFP/OD defaults; every instrument-specific value is flagged.

Note: this DNA-free copper screen reads **OD600 only** (no fluorophore is present without DNA). The sfGFP settings are specified now so the same plate protocol carries directly into Wet EXP-1/M1 once the PCu-B0032-sfGFP construct is in hand.

| setting | proposed value | provenance |
|---|---|---|
| OD600 | 600 nm absorbance, path-length corrected if available | standard |
| sfGFP excitation | 485 nm (or 470-490 filter) | vendor-manual default for sfGFP/GFPmut |
| sfGFP emission | 510-515 nm | vendor-manual default |
| gain | fixed manual gain, set once on a bright control and held across the plate | **requires Prof. Nishida confirmation** (instrument-specific) |
| read mode | endpoint + optional kinetic (OD every 15-30 min for growth curves) | standard |
| shaking | double-orbital, 3-5 s before each read | standard |
| temperature | 37 C | standard |
| blank wells | column 11 (medium-only), subtract per medium | this layout |
| background subtraction | subtract medium-only blank from every OD and every GFP read | standard |

Label on the final settings sheet: vendor-manual values (ex/em, read mode) vs values still needing Prof. Nishida's sign-off (gain strategy, exact filter set for the installed instrument).

---

## 3. Replicates and 96-well plate layout

- **Biological replicates:** 3 independent overnight cultures of untransformed DH5alpha per ladder.
- **Technical replicates:** each concentration has a single well per biological replicate; the instrument's internal averaging is one well read once, not a technical replicate. A bad well (bubble, edge drift, pipetting error) corrupts that point with no within-plate backup, and an anomalous point cannot be adjudicated within the plate. For a range-finding screen, n=3 biological replicates is defensible; if bench capacity allows, add a second technical well per concentration (drops the plate to 2 ladders with fewer points, or 1 ladder with duplicate wells).
- **Controls:** column 11 = medium-only blank (no cells); column 12 = 0-Cu growth control (cells, no copper). The ladder's own first point (0 Cu) is also a 0-Cu growth control, so col 12 is a deliberate second baseline; repurpose it to a within-plate positive-inhibition control (a known-inhibitory dose) if a second baseline is not wanted.
- **Layout (see figure):**
  - Rows A-C: LB ladder x 3 biological reps (cols 1-10), blank (11), 0-Cu control (12).
  - Rows D-F: M9 ladder x 3 biological reps (same column scheme).
  - Rows G-H: medium-only evaporation guard.
- **Edge effects.** The G-H guard covers only the bottom edge; row A and the outer columns (col 1 = 0-point, col 12 = control) still sit on the plate perimeter, the wells most prone to evaporation and edge drift. A full-perimeter guard (fill row A, row H, col 1, col 12 with medium) conflicts with using all 12 columns for data plus controls. Two defensible options: (i) accept the trade-off and **rely on within-block randomisation** of well-to-culture assignment to decouple position from concentration, or (ii) drop one biological replicate to free a perimeter row/column for a full guard. State the choice on the bench sheet.

---

## 4. Model linkage (explicit)

- The DNA-free copper screen **constrains the safe/usable copper concentration range** for Wet EXP-1 and the future Wet M3 copper recorder: the recorder must operate below the host growth-inhibition onset, or growth confounds the readout.
- It is **not** a measurement of Bxb1 flipping, integrase leak, or pBAD leakage. Those require the DNA constructs (Wet EXP-2) and are out of scope for the DNA-free week.
- Output feeds back into the models as an **upper bound** on the copper doses the recorder can be characterised at, and as a host-fitness anchor for the Module 4 fitness-cost parameter c.

---

## Deliverable checklist (Section 4)
- [x] CuSO4 growth-tolerance ladder, justified, not the recorder ladder (Sec 1)
- [x] Plate-reader settings, vendor-vs-confirmation labelled, instrument model flagged TBD (Sec 2)
- [x] Replicates + 96-well layout with blanks and untransformed controls (Sec 3)
- [x] Model-linkage statement, with the "not a flip/leak measurement" caveat (Sec 4)

---

## Sources
- Franke S, Grass G, Rensing C, Nies DH. 2003. J Bacteriol 185:3804-3812.
- Rosenberg M, Umerov S, Teaer CM, Ivask A. 2025/2026. Microbiol Spectrum (biorxiv 2025.08.27.672559).
- Macomber L, Imlay JA. 2009. PNAS 106:8344-8349.
- Pang Y et al. 2020. Front Microbiol 10:3031 (WMC-007 copper bioreporter).
- M2 v3b Revision Instructions, Section 9 (WMC-007 strain background).
