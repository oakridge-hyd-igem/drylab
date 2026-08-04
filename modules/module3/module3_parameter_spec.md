# Module 3 — Stochastic Flip Model: Parameter Specification

**BioChronos** · Dry lab / modelling · *specification only — no simulation code yet*

This document defines the model on paper. Its parameter table is the direct input to
**Module 7** (one-at-a-time sensitivity sweep): the "Plausible range" column **is** the
Module 7 sweep range. All parameters are collected here so that a single `config` dict
(code convention) can be passed into the simulation and wrapped by the sweep without
refactoring.

**Core question the model must answer:** *what is the minimum pulse duration that still
flips a detectable (≥5%) fraction of cells?* Every parameter below is chosen with that
question in mind.

---

## 1. Model specification

### 1.1 What we track and why

Under the project's QSSA convention (mRNA at steady state, `d[mRNA]/dt = 0`, mRNA
turns over far faster than protein) we do **not** model mRNA explicitly. We model
integrase **protein** accumulation directly, and the flip stochastically. The model
tracks three processes:

1. **Integrase production** — driven by the sensor promoter, a Hill function of analyte
   concentration `A`, plus a basal *leak* term at `A = 0`.
2. **Integrase degradation** — dilution by cell growth (always present) **plus** active
   degradation from an LVA ssrA tag where the design includes one.
3. **The flip** — unflipped DNA → flipped DNA. Irreversible, one-way, per cell,
   catalysed by integrase. This is the recorded event.

The flip **must** be stochastic. It is discrete and irreversible in each cell, and the
readout is the **population flipped fraction**, so a deterministic ODE for "fraction
flipped" would hide exactly the cell-to-cell variability that decides whether a short
pulse crosses the detection threshold. Production and degradation are smooth and
high-copy, so we treat them deterministically per cell; only the flip is drawn
stochastically. This "deterministic protein trajectory + stochastic flip hazard" hybrid
is the simplest model that answers the core question — it is an agent-based population of
`N` cells, each carrying its own integrase level and flip status.

### 1.2 Per-cell state variables

| Variable | Meaning | Type |
|---|---|---|
| `I(t)` | integrase protein count in the cell | non-negative real (copies/cell) |
| `s(t)` | flip status: `0` = unflipped, `1` = flipped | binary, latches at 1 (irreversible) |

Population state is the set `{ (I_i, s_i) : i = 1..N }`. The readout is the flipped
fraction `F(t) = (1/N) Σ s_i`.

### 1.3 Governing equations

**Analyte input (the pulse):**

```
A(t) = A_pulse   for 0 ≤ t ≤ t_pulse
A(t) = 0         otherwise
```

**Integrase production (per cell), Hill activation by analyte:**

```
dI/dt = β(A)  −  γ · I

β(A) = β_leak  +  β_max · A^n / (K^n + A^n)
```

`β(A)` is the QSSA-collapsed production rate (transcription × translation folded into one
protein-production term). `β_leak` is the basal rate at `A = 0` — the false-positive
driver, and the term the leak-suppression strategies (weak RBS, LVA tag) act on.

**Integrase degradation:**

```
γ = γ_dil            (untagged: dilution by growth only)
γ = γ_dil + γ_LVA    (LVA-tagged: dilution + active degradation)
```

**The flip (stochastic, per unflipped cell):** each cell with `s = 0` flips as a Poisson
process with an integrase-dependent hazard `λ(I)`:

```
λ(I) = k_flip · I / (I + K_I)          (saturating in integrase)

P(flip during [t, t+dt] | s=0) = 1 − exp(−λ(I(t))·dt) ≈ λ(I(t))·dt
```

Once `s` flips to 1 it never returns (irreversibility). The saturating (Michaelis–Menten)
form encodes that the flip needs a synaptic complex (a Bxb1 tetramer bridging attB/attP):
below `K_I` integrase the flip rate rises with integrase; above it, the rate plateaus at
`k_flip`.

**Readout:** at readout time the event is called *recorded* if `F ≥ F_detect` (default
0.05, flow cytometry).

> **Simulation loop (for Module 4, stated here so the parameters map cleanly):** for each
> timestep `dt`, (i) advance every cell's `I` by the production/degradation equation,
> (ii) for each unflipped cell draw a uniform random number and flip it if
> `u < 1 − exp(−λ(I)·dt)`, (iii) record `F(t)`. Requires `dt ≪ 1/γ` and `dt ≪ 1/k_flip`.

---

## 2. Parameter table

Values in **per-hour** time units and **protein copies per cell** throughout (state the
unit once, keep it consistent). "ESTIMATED" marks any value with no direct measurement —
Module 7 exists precisely to test how much each of these actually matters, so we give an
order-of-magnitude value and a wide range rather than stalling on precision.

> **Baseline provenance (v5):** the flip and clearance baselines (`k_flip`, `γ_dil`,
> `K_I`, `β_leak`) are anchored to the fitted parameters of Hsiao et al. 2016 — the one
> published Bxb1 *recording* model that fit population data — rather than to textbook
> order-of-magnitude defaults. This makes the baseline result literature-matched and
> conservative. The Alon textbook values (γ = ln2/τ ≈ 1.4 h⁻¹, leak ~10⁻³) remain inside
> the sweep ranges, so Module 7 still spans them.

### Hazard-form sensitivity (validation against Hsiao's mechanistic model)

The flip hazard above is **first-order saturating**, λ(I) = k_flip·I/(I+K_I). This is a
simplification of Hsiao 2016's mechanistic form. Hsiao Eq. (2) models the flip as a
**4th-order tetramer** binding function — a flip requires four integrase monomers to
occupy attB/attP and form the synaptic tetramer, so λ(I) = 0 for I<4 and rises steeply
thereafter. The two forms were compared directly (see `module3_validation.md`, which
implements Hsiao's exact equation from the fetched full text):

- **They agree on genuine events.** At induced (high) integrase both hazards give the
  same near-complete flip, so the core claim — a real pulse is recorded, and the
  recorder is leak-limited not sensitivity-limited — holds under both and does not
  depend on the hazard form.
- **They diverge only at low integrase, where the first-order form is CONSERVATIVE.**
  The tetramer needs ≥4 integrase, so at the leak steady-state levels of the
  suppressed designs the first-order model *overestimates* false positives (e.g.
  weak RBS: 91% first-order vs ~0% tetramer at 24 h; LVA tag: 98% vs 35%). For a
  false-positive-averse recorder this is the right (pessimistic) default to headline —
  the suppressed designs are only *safer* under the more mechanistic tetramer.

**Status:** the first-order hazard is retained as the conservative baseline for all
Module 3 headline results. The Hsiao tetramer is a named **structural-sensitivity
variant for Module 7** (a model-form axis, alongside the numeric parameter sweeps).
Both are valid datapoints that bracket the flip-kinetics assumption; neither supersedes
the other. Full detail and figures: `module3_validation.md`.


| Symbol | Meaning | Unit | Baseline value | Plausible range (Module 7 sweep) | Source |
|---|---|---|---|---|---|
| `β_max` | Max integrase production rate (promoter fully induced) | copies·cell⁻¹·h⁻¹ | 1000 | 100 – 10 000 | Alon *Intro to Systems Biology*, Ch. 2 (strong-promoter steady state ~10³–10⁴ copies ⇒ β_max ≈ γ·I*); ESTIMATED magnitude |
| `n` | Hill coefficient of sensor promoter | dimensionless | 2 | 1 – 4 | Alon Ch. 2 (cooperative TF activation typically n = 1–4) |
| `K` | Half-max (activation) analyte concentration | µM | 5 | 1 – 15 | **Analyte LOCKED = copper.** CueR/P(copA) sensor; K grounded in Module 1 (copper-biosensor operating windows, Pang 2020 / Stoyanov 2001 abstract) |
| `β_leak` | Basal integrase production at A = 0 (leak) | copies·cell⁻¹·h⁻¹ | **100 (= 10%·β_max, T1-grounded to copper)** | 10 – 500 (1% floor – 50% worst natural) | T1: natural CueR/P(copA) fold-change 7–18× (Fu 2024, full text). Hsiao's 1–2% integrase-construct leak is the sweep floor. Old baseline was 10 (=1%); see DECISIONS.md / T1_resolution_beta_leak.md. Original Hsiao fit leak = 1–2% of production. (Alon Ch. 2 textbook default ~10⁻³·β_max is the low end of the range.) **leak-suppression target** |
| `γ_dil` | First-order integrase clearance, untagged (dilution + slow degradation) | h⁻¹ | **0.3 (t½ ≈ 2.3 h)** | 0.3 – 2.1 | Hsiao et al. 2016 fitted k_deg = 0.3 h⁻¹ (effective clearance in their growth conditions). Upper end 2.1 = fast dilution at τ = 20 min (Alon, γ = ln2/τ) |
| `γ_LVA` | **Extra** degradation from LVA ssrA tag (added to γ_dil) | h⁻¹ | 1.04 (t½ ≈ 40 min) | 0.7 – 1.7 (t½ ≈ 60–25 min) | Andersen et al. 1998 (GFP-LVA t½ ≈ 40 min in *E. coli*). Untagged γ = γ_dil = 0.3 h⁻¹; **LVA-tagged γ = γ_dil + γ_LVA ≈ 1.34 h⁻¹** (effective t½ ≈ 31 min) |
| `k_flip` | Max per-cell flip hazard (integrase saturating) | h⁻¹ | **0.4** | 0.1 – 10 | Hsiao et al. 2016 fitted k_flip = 0.4 h⁻¹ (single flip term, no saturation). Zhang et al. 2022 show Bxb1 rate tunable over orders of magnitude by attP → wide range |
| `K_I` | Integrase copies at half-max flip rate | copies/cell | **10** | 1 – 1000 | Hsiao et al. 2016 fitted Kd = 10 molecules (tetramer/synaptic-complex requirement; Bxb1 acts as a dimer-of-dimers); ESTIMATED upper range |
| `A_pulse` | Analyte concentration during pulse (**input**) | µM | 10 (≈ 10·K, saturating) | 0 – 100 | Input variable; swept to build dose–duration surface |
| `t_pulse` | Pulse duration (**input — the quantity we solve for**) | h | scanned | 0.05 – 24 | Input variable; the primary sweep axis |
| `F_detect` | Detection threshold (flipped fraction called positive) | fraction | 0.05 | 0.01 – 0.50 | Project convention (flow cytometry ~5%); 0.01 = sensitive flow, ≥0.1 = bulk PCR orientation assay |
| `N` | Population size (agents simulated) | cells | 2000 | 100 – 100 000 | Simulation parameter; N sets the finest resolvable F (≥ 1/N) and the sampling noise on F (~√(F(1−F)/N)) |
| `dt` | Simulation timestep | h | 0.01 (≈ 0.6 min) | 0.001 – 0.05 | Numerical; must satisfy dt ≪ 1/γ and dt ≪ 1/k_flip |

**Derived quantities (not swept independently — computed from the above):**

- Steady-state integrase at full induction: `I* = (β_leak + β_max)/γ`. Untagged ≈ 3370
  copies (γ = 0.3); LVA-tagged ≈ 754 copies (γ = 1.34) — baseline values.
- Steady-state leak integrase: `I_leak = β_leak/γ` ≈ **333 (untagged) / 74.6 (LVA)** copies
  at the grounded β_leak=100 (was 33 / 7.5 at the old 1% leak) — the false-positive
  pressure. Note LVA still cuts steady-state leak integrase ~4×. **Grounding the leak 10×
  shortens the trust horizon** (time before leak alone crosses the 5% detection threshold):
  untagged 0.61→**0.26 h**, LVA 0.67→**0.27 h**, weak RBS 1.84→**0.61 h**, combined
  2.59→**0.67 h** — see module3_results.csv (regrounded).
- Integrase response/decay time constant: `1/γ` ≈ 200 min (untagged) / 45 min (LVA-tagged).

---

## 3. Flagged uncertainties

Ordered by how much they need Module 7's attention. The first three are both the least
certain **and** sit directly on the false-positive / false-negative tension, so they
matter most.

1. **`k_flip` — flip rate constant (LEAST confident, highest priority).** No clean in-vivo
   per-cell hazard exists in the literature; measured Bxb1 kinetics are bulk in-vitro or
   whole-population qPCR/fluorescence, and Zhang et al. (2022) show it is tunable across
   orders of magnitude by attP sequence alone. This single parameter sets the minimum
   pulse duration almost linearly — a 10× error here is a 10× error in the answer. Sweep
   the full 0.1–10 h⁻¹ range; expect it to dominate the sensitivity ranking.

2. **`β_leak` — leak / basal production (drives false positives).** Basal promoter
   activity is notoriously construct- and context-dependent; our 10⁻³·β_max is a
   textbook default, not a measurement. It is the entire reason weak-RBS / LVA / combined
   leak-suppression strategies exist. Critically, it is in **tension** with `k_flip` and
   `γ`: suppressing leak to kill false positives also blunts the response to short pulses.
   Module 7 should sweep `β_leak` jointly against `t_pulse` to map that trade-off.

3. **`K_I` — integrase level for half-max flip (shape of the flip's integrase
   dependence).** ESTIMATED. Controls whether a brief, low integrase transient is enough
   to flip. If `K_I` is high, short pulses fail even when leak is low — directly a
   false-negative lever. Interacts with `k_flip` and `γ`.

4. **`γ_LVA` and the untagged-vs-tagged `γ` gap.** The LVA half-life (~40 min) is
   well-sourced (Andersen 1998) but tag efficiency varies with ClpXP load and growth
   phase; faster tags (LAA, ~minutes) sit outside our range and would be a separate
   design point. `γ` controls both how fast integrase clears after the pulse ends (memory
   fidelity) and steady-state leak integrase, so it appears on both failure modes.

5. **`K` and `n` — sensor promoter parameters.** The analyte is **LOCKED = copper** (CueR/P(copA)); `K` is grounded in Module 1 (copper-biosensor operating windows, baseline 5 µM, range 1–15 µM) and `n≈1` is the grounded sensor Hill coefficient (Stoyanov 2001 abstract + MerR single-site mechanism; T2 resolved — see Module 1). Module 3's flip results are n-independent for supra-threshold events. These shift the dose axis but, because `A_pulse` is
   set saturating (≈10·K) for the minimum-duration question, they matter less for the core
   result than for a dose–response study. Revisit once the sensor is chosen.

Lower concern (standard textbook values, or numerical): `β_max`, `γ_dil`, `F_detect`, `N`,
`dt`. `β_max` mostly rescales `I` and is partly absorbed by `k_flip`/`K_I`; `N` and `dt`
are convergence controls, not biology.

---

## References

- U. Alon, *An Introduction to Systems Biology: Design Principles of Biological Circuits*
  (Chapman & Hall/CRC), Ch. 1–2 — baseline *E. coli* production, dilution (γ = ln2/τ), and
  Hill-function parameters.
- J. B. Andersen et al. (1998) "New unstable variants of green fluorescent protein for
  studies of transient gene expression in bacteria," *Appl. Environ. Microbiol.* 64:2240–2246
  — ssrA LVA/AAV/ASV degradation-tag half-lives (LVA ≈ 40 min).
- Hsiao, Hori, Rothemund & Murray (2016) "A population-based temporal logic gate for
  timing and recording chemical events," *Mol. Syst. Biol.* 12(5), 2016,
  doi:10.15252/msb.20156663 — Markov model of integrase flipping in a bacterial
  population; step/pulse inputs. **Directly analogous recording model.** Fitted parameters
  (from the paper): k_flip = 0.4 h⁻¹, k_deg = 0.3 h⁻¹ (t½ ≈ 2.3 h), Kd = 10 molecules,
  k_leak = 1–2% of k_prod; no cooperativity observed in Bxb1 DNA binding. These are the
  primary cross-check values for our `k_flip`, `γ`, `K_I`, and `β_leak`.
- Zhang, Azarin, Sarkar et al. (2022) "Model-guided engineering of DNA sequences with
  predictable site-specific recombination rates," *Nat. Commun.* 13:4152,
  doi:10.1038/s41467-022-31538-3 — qPCR-measured Bxb1 inversion rate, tunable over orders
  of magnitude by attP sequence. (Basis for the wide `k_flip` sweep range.)
- Friedland, Lu, Wang et al. (2009) *Science* 324:1199–1202; Bonnet, Subsoontorn & Endy
  (2012) "Rewritable digital data storage in live cells via engineered control of
  recombination directionality," *PNAS* 109:8884–8889; Siuti, Yazbek & Lu (2013)
  *Nat. Biotechnol.* 31:448–452 — integrase-based counting/logic and irreversible memory
  (context for the recording concept).
- Chao, Travis & Church (2021) "Measurement of large serine integrase enzymatic
  characteristics in HEK293 cells," *FEBS J.* 288:6410–6427, doi:10.1111/febs.16037 —
  large-serine-integrase catalytic rate and substrate affinity.
