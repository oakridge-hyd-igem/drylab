# Module 5: Multi-State Separability

**BioChronos (Oakridge-HYD iGEM 2026) · analyte: COPPER · 2026-07-27**

Maps to the 0 / 1 / 2 multi-state design question (plan Decision 0, C2, D7).

## 1. Question

The recorder should ideally report not just *whether* a copper event happened
(binary) but *how much*, a low/transient exposure (**State 1**) vs a high/sustained
one (**State 2**), on top of no exposure (**State 0**). A single cell either flips or
not, so at the population level "State 1" means a **partially flipped population**.
The plan defines a population-fraction threshold (≈20–60% flipped = State 1, >60% =
State 2) and asks: **does single-integrase threshold tuning give separable states, or
is a two-integrase architecture required?** Metric: distribution overlap / 3-state
classification accuracy.

## 2. Model (built on the Module 3 engine)

A population readout (flow cytometry, or gel/qPCR of the flipped locus) reports a
flipped **fraction**. Across replicate populations that fraction spreads from two
grounded sources:

1. **Extrinsic expression noise**: population-to-population variation in integrase
   expression, modelled as a log-normal scaling of β with **CV = 0.30**. This is the
   measured *E. coli* extrinsic-noise floor: Taniguchi et al. 2010 report the
   high-expression protein-noise plateau at η_p² ≳ 0.1 ⟹ CV ≈ √0.1 ≈ 0.3. **This
   single parameter controls the whole verdict, and it is grounded, not assumed.**
2. **Finite-sample binomial noise**: N cells sampled per population (default 2000).

For each of the three canonical exposures we simulate R=1500 replicate populations,
build the flipped-fraction distribution, and score an optimal two-threshold 3-state
classifier. Flip kinetics (k_flip=0.4/h, K_I=10), leak (β_leak grounded to 10% of β_max, T1),
and the copper sensor (K=5 µM, n=1) are inherited from Modules 1/3. The combined
design's induced level β_max=100 follows **Module 1's convention that a weak RBS
(×0.1) scales induced AND basal expression equally** (β_max = BETA_MAX×0.1 = 100,
β_leak = 100×0.1 = 10); the untagged design keeps β_max=1000, β_leak=100. Module 5
adds the population noise model (extrinsic CV + binomial sampling); it introduces no
new flip-kinetics parameters.

Canonical exposures: **State 0** = none (leak only); **State 1** = 2 µM for 0.5 h
(low/transient); **State 2** = 20 µM sustained through the readout.

## 3. Results

### 3.1 A single integrase separates three states: if leak is suppressed

At the grounded extrinsic noise (CV=0.3) and the optimal readout (**1.25 h**), the
leak-suppressed **combined** design resolves all three states at **94% accuracy**,
with small State-1/State-2 distribution overlap (0.06):

| Readout | mean S0 | mean S1 | mean S2 | overlap(S1,S2) | 3-state accuracy |
|---|---|---|---|---|---|
| 0.5 h | 0.03 | 0.08 | 0.11 | 0.19 | 91.8% |
| 1.0 h | 0.09 | 0.17 | 0.24 | 0.11 | 93.8% |
| **1.25 h** | **0.12** | **0.21** | **0.30** | **0.06** | **94.1%** |
| 1.5 h | 0.15 | 0.25 | 0.36 | 0.05 | 93.2% |
| 2.0 h | 0.22 | 0.31 | 0.46 | 0.02 | 90.9% |
| 3.0 h | 0.33 | 0.42 | 0.62 | 0.01 | 86.0% |

There is an **optimal readout window (~1–1.5 h)**: too early and the states have not
yet spread apart; too late and leak drags State 0 up into State 1, eroding accuracy.

### 3.2 Leak suppression is the main lever; a second integrase adds a modest gain

The **untagged** design (grounded 10% leak) reaches only **79%** three-state accuracy
at the same readout: leak floods the population, collapsing State 0/1/2 toward each
other (means 0.32 / 0.36 / 0.38, overlap 0.49). Suppressing leak (weak RBS + LVA tag)
is what lifts accuracy to 94% (+15 points), and it is exactly the design Modules 1/2/4
already select for a different reason.

**On the two-integrase alternative (D7/D8):** a second orthogonal integrase gives a
2-bit genotype (locus A behind a low-threshold sensor, locus B behind a high-threshold
one) rather than a point on a continuous fraction axis. Simulated properly (both loci
on the combined chassis, K_A=2 µM / K_B=15 µM, optimal 2D classifier), it reaches
**~98%** at the same 1.25 h readout and grounded noise, i.e. **~4 points above the
single-integrase design**. So the second integrase does give a real, if modest,
separability gain. The trade-off is construct complexity, orthogonality validation
(low Bxb1 cross-reactivity, D8) and a second att-site design.

**Verdict:** the single leak-suppressed integrase already clears a usable 3-state
readout (~94%); the second integrase buys ~4 points at meaningful added complexity. For
the *binary* recorder (the project's core claim) a second integrase is **not required**;
for a high-confidence *multi-state* readout it is a worthwhile but optional upgrade whose
value should be weighed against the measured expression noise (§3.3) and build cost.

### 3.3 The verdict depends on expression noise (grounded, but flagged)

| Extrinsic CV | 0.1 | 0.2 | **0.3 (grounded)** | 0.5 | 0.75 | 1.0 |
|---|---|---|---|---|---|---|
| 3-state accuracy (combined, 1.25 h) | 99.7% | 98.4% | **93.3%** | 80.5% | 68.5% | 60.9% |

(The CV=0.3 point here reads 93.3%, and §3.1's readout sweep reads 94.1% at the same
1.25 h / CV=0.3 condition. These are two independent Monte-Carlo runs of the *same*
condition, differing in random seed (readout sweep seeds 101-103, R=1500; noise sweep
seeds 201-203, R=1200). Re-running confirms the gap is sampling scatter: at a fixed seed
set the two population sizes agree to ~0.2 points, and the residual ~0.8-point spread is
driven by the independent seeds, not by any real effect. Both round to "usable ~94%.")

Separability is strong at and below the measured *E. coli* noise floor (CV≈0.3) but
degrades if the actual construct shows larger population-to-population variability.
If wet-lab flow cytometry reveals CV > ~0.5, the single-integrase 3-state readout
would need either a larger sampled N, a tighter readout window, or a move to the
two-integrase genotype architecture (which retains its ~4-point margin). **This is the
one thing to check in EXP3.**

## 4. Recommendation to the wet lab

- **Build the leak-suppressed combined construct** (single integrase): it is the same
  construct Modules 1/2/4 recommend, and it delivers ~94% 3-state separability at
  grounded noise. This suffices for the binary recorder and gives a usable multi-state
  readout.
- **A second orthogonal integrase is an optional upgrade**, not a requirement: it adds
  ~4 points (~98%) for high-confidence multi-state discrimination, at the cost of a
  second att-site/integrase and orthogonality validation (D7/D8). Decide it against the
  measured expression noise and build budget, not the separability requirement alone.
- **Read out at ~1.25 h** after the exposure for best State-1/State-2 discrimination.
- **In EXP3, measure the population-to-population CV of the flipped fraction** (across
  replicate cultures / flow-cytometry runs). That number confirms or refutes the
  separability verdict; our CV=0.3 is a literature floor, not a device measurement.

## 5. Caveats

1. **Separability rests on extrinsic noise CV=0.3**: grounded to the *E. coli*
   noise floor (Taniguchi 2010), but the real construct could be noisier. Swept in §3.3
   and routed to Module 7.
2. **Deterministic per-cell integrase trajectory**: the single-cell flip probability
   uses the shared Module 3 mean trajectory; intrinsic single-cell integrase
   fluctuation is folded into the extrinsic-CV term rather than modelled separately.
3. **State definitions are exposures, not ground truth**: "low/transient" (2 µM,
   0.5 h) and "high/sustained" (20 µM) are representative; the real State 1/2 copper
   levels are a design choice tied to the deployment target.
4. **First-order flip hazard** (conservative vs Hsiao tetramer, per Module 3
   validation), the tetramer form would sharpen the high/low contrast, if anything
   improving separability.

## 6. References

- Taniguchi, Choi, Li, Chen, Babu, Hearn, Emili & Xie (2010) "Quantifying *E. coli*
  proteome and transcriptome with single-molecule sensitivity in single cells,"
  *Science* **329**:533–538, doi:10.1126/science.1188308, extrinsic protein-noise
  floor η_p² ≳ 0.1 (CV ≈ 0.3). CrossRef + full text (PMC) verified this session.
- Elowitz, Levine, Siggia & Swain (2002) "Stochastic gene expression in a single
  cell," *Science* **297**:1183–1186, doi:10.1126/science.1070919, intrinsic vs
  extrinsic noise decomposition. CrossRef-verified.
- Flip kinetics, leak and sensor parameters inherited from Modules 1 & 3 (Hsiao 2016;
  Fu 2024 copper-promoter grounding; copper sensor CueR/P(copA)).
