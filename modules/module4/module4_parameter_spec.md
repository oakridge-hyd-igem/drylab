# Module 4 — Irreversibility / Memory Stability

**BioChronos** · Dry lab / modelling · Oakridge-HYD iGEM 2026

**Maps to:** the memory-persistence claim — a written state must survive after
the signal disappears — and the EXP2 Memory Stability wet-lab test.

**Core question:** once the Bxb1 flip has recorded a state, does that state persist
through continued growth, cell division, and dilution over a readout delay, or
does it decay?

---

## 1. Why DNA memory is different

The record is a **DNA inversion**, so it is **heritable** — when a flipped cell
divides, both daughters inherit the flipped chromosome. This is the fundamental
contrast with a protein/metabolite memory: once production stops, a protein
memory **halves every generation** by dilution, whereas the DNA flip does not
dilute at all. This is the modelling justification for choosing an integrase
recorder over a protein-based one.

The only two processes that can erode a DNA record are:
1. **Fitness cost** `c` — if the flipped state slows growth, unflipped cells
   outcompete flipped ones and the population fraction drifts down.
2. **Continued leak** `λ` — but because the flip is one-way, leak *increases* the
   flipped fraction; it can never erase it. Leak is a false-positive problem
   (Module 2), not a memory-loss problem.

## 2. Model

Two subpopulations, U (unflipped) and F (flipped, irreversible), with growth rate
`μ` set equal to the dilution rate `γ_dil = 0.3 h⁻¹` from Modules 2/3
(doubling time 2.31 h):

```
dU/dt = μ·U − λ·U
dF/dt = μ·(1−c)·F + λ·U
```

The flipped fraction `f = F/(U+F)` has a closed form (validated to machine
precision, ≤7×10⁻⁸, against the full two-population ODE):

```
df/dt = (1−f)·(λ − μ·c·f)
```

- **Equilibrium:** `f* = min(1, λ/(μ·c))` — approaches 1 as cost → 0.
- **Cost-only decay** (λ=0): the odds `f/(1−f)` decay as `2^(−c·g)` over `g`
  generations, so the **memory odds half-life is 1/c generations**.

Leak hazards `λ` per design are inherited from Module 2, at the T1-grounded
β_leak=100 (10% copper-promoter leak): untagged 0.388, weak RBS 0.308, LVA tag
0.353, combined 0.171 h⁻¹. (At the old 1% leak these were 0.308/0.100/0.171/0.028;
grounding raises every hazard, which only reinforces the flipped record.)

## 3. Results

### 3.1 The memory is robust — it persists and even self-reinforces

Starting from a genuine recording event (~90% flipped, per Module 3) and holding
for 48 h under the combined (lowest-leak) design:

*(β_leak grounded to 10%, T1 — combined-design leak hazard λ=0.171 h⁻¹.)*

| Fitness cost | f at 0 h | f at 24 h | f at 48 h |
|---|---|---|---|
| 0% | 0.90 | 1.00 | 1.00 |
| 2% | 0.90 | 1.00 | 1.00 |
| 5% | 0.90 | 1.00 | 1.00 |
| 10% | 0.90 | 1.00 | 1.00 |

The flipped fraction **rises to ~100% within a day** and never falls, for every
fitness cost up to 10%. At the grounded (10×-higher) copper leak, residual leak
reinforces the record far faster than even a large fitness cost can erode it —
more strongly than at the old 1% leak, where a 10% cost held the fraction flat at
0.90 rather than driving it up. Memory is emphatically not the binding constraint.

### 3.2 Memory half-life is long compared to a readout

For the worst realistic case (cost-only decay, no reinforcing leak):

| Fitness cost | odds half-life (generations) | odds half-life (hours) |
|---|---|---|
| 2% | 50 | 116 h |
| 5% | 20 | 46 h |
| 10% | 10 | 23 h |

Even a **substantial 10% growth penalty gives a ~23 h memory odds half-life** —
long relative to a same-day or next-day readout. A more realistic ≤5% cost gives
≥46 h. The record is safe across the EXP2 stability-test timeframe.

### 3.3 Long-run equilibrium

`f* = min(1, λ/(μc))`. At the T1-grounded leak (β_leak=100), **every design —
including the leak-suppressed combined design — reaches `f* = 100%`** across the
entire plausible cost range: the smallest hazard is the combined design's
λ=0.171 h⁻¹, and `λ/(μc) ≥ 1` holds until the fitness cost exceeds **c ≈ 57%**
(`c > λ/μ = 0.171/0.3`), far beyond any realistic value. So at grounded copper
leak there is **no memory/false-positive trade-off left** — the leak-suppression
that helps false positives (Module 2) no longer measurably weakens memory
self-reinforcement, because even the suppressed leak is strong enough to pin the
record at 100%. (At the old 1% leak, the combined hazard was λ=0.028 and this
equilibrium did dip to f*=92.6% at c=10%; grounding the leak 10× higher erased
that residual trade-off.)

## 4. The Module 2 ↔ Module 4 trade-off (key cross-module result)

Leak plays **opposite roles** in the two failure modes:
- **Module 2 (false positive):** leak is the enemy — it flips cells with no event.
- **Module 4 (memory):** leak is a friend — it reinforces a genuine flip against
  fitness-cost erosion.

The combined leak-suppressed design is correctly optimised for **Module 2**
(false positives dominate the risk for an irreversible recorder), and at the
T1-grounded copper leak it pays **no** memory-stability price: it rises to ~100%
retained by 48 h even at 5–10% fitness cost (its residual leak still reinforces
the record faster than cost erodes it). Memory stability is **not** a binding
constraint on the design; false-positive control is.

## 3.4 Stochastic DNA model and validation against measured data

The continuous fraction model above assumes an idealised single stable locus. To
test the record against real failure modes — finite-population **genetic drift**
and inheritance through division — the flip state is also modelled as a
DNA-borne "allele" in a **Wright-Fisher population** (`wright_fisher()` in
`module4_memory_model.py`): each generation applies selection (flipped cost `c`),
one-way leak conversion (U→F; the flip is irreversible so there is no reverse
mutation), and binomial resampling of `N` offspring (drift).

**Fitness cost — estimate informed by the literature (not a derived value).** Two
measurements inform the plausible magnitude: (i) increased Bxb1 expression imposed
**no measurable metabolic burden** on the control construct (Gonzalez-Colell &
Macia 2025, *Sci. Rep.* **15**:31311, full text) — i.e. integrase burden ≈ 0; and
(ii) the expression–growth tradeoff is approximately linear, and a *strong*
constitutive reporter (eGFP fraction ≥0.07 mg/mg DCW) is needed to push growth
down to 0.59 h⁻¹ (Bienick et al. 2014, *PLoS ONE* **9**:e109105, full text) —
a single moderate reporter sits well below that. **These bound the cost as
*small*, but neither paper yields a single number for this construct.** The values
used here — **nominal baseline c = 2%, sweep range 0–5%** — are therefore a
*qualitative-informed estimate*, not a quantity derived from either paper. They
replace the earlier arbitrary 0–10% sweep with a range consistent with
"burden is small," and the module's conclusion is checked for robustness across
the whole range below (it does not hinge on the exact value).

**Robustness to the cost estimate.** Because c is an estimate, the conclusion is
checked across the whole plausible range (combined design, residual leak on,
flipped fraction at 90 generations): c=0% → 100.0%, c=2% → 100.0%, c=5% → 100.0%. At the grounded 10% copper leak the reinforcement dominates any plausible fitness cost, so retention pins at ~100% across the entire 0–5% range — even more robustly than at the old 1% leak (which gave 97.5–100%). So the
"record survives over many generations" conclusion does not depend on the specific
value chosen.

**Comparison to a literature benchmark (caveated).** Integrase/recombinase memory
is widely reported to persist over many cell generations (Siuti, Yazbek & Lu 2013,
*Nat. Biotechnol.* 31:448–452). *Caveat: the specific figure often quoted (~"90+
generations") could NOT be confirmed from the paper's full text or abstract this
session — the article is closed-access and no abstract was retrievable. It is
recorded here as a literature-reported order-of-magnitude target, not a
session-verified measurement.* The model, run to 90 generations at a nominal 2%
cost, gives:

| Design | flipped @ 90 gen (model, c=2%) |
|---|---|
| untagged | 100.0% |
| combined | 99.8% |

i.e. the model predicts retention over the many-generation timescale that
recombinase recorders are reported to achieve — a consistency check, not a
quantitative fit to a verified datapoint. Record loss appears
**only** in an artificial corner — leak fully suppressed (λ=0) *and* high cost
(≥5%) *and* small founding population — where drift can drive the flipped allele
toward loss (e.g. leak-free, c=5%: 8% retained at 90 gen). Under realistic
conditions (any residual leak, grounded cost) retention is >97% at 90 generations
across all designs.

## 5. Caveats

1. **Fitness cost c is a literature-informed estimate, not a derived or measured
   value.** Nominal baseline 2%, sweep 0–5%. The burden papers (Gonzalez-Colell &
   Macia 2025; Bienick 2014) establish only that the cost is *small*; they do not
   pin a number for this construct. The conclusion (memory retained over many
   generations) holds across the entire 0–5% range — see robustness check — so it
   does not depend on the specific value chosen. It is NOT a quantitative fit to a
   verified benchmark.
2. **Deterministic population fractions.** Segregation/partitioning noise at
   division is not modelled; at population scale (Module 3 uses N=2000) this is
   negligible, but a small founding population could drift.
3. **No plasmid loss / mutational inactivation.** Assumes the flipped locus is
   chromosomal and genetically stable over the timescales modelled (hours–days),
   consistent with a genomically integrated device.
4. **Growth rate = dilution rate assumption.** μ is set to γ_dil for consistency;
   in a real deployment (stationary phase, encapsulation) growth may be much
   slower, which only *strengthens* memory (less cost-driven washout). Swept in
   Module 7.

## 6. References

- U. Alon, *An Introduction to Systems Biology* — growth/dilution, two-species
  competition dynamics.
- Bxb1 irreversibility and heritable state: Bonnet, Subsoontorn & Endy (2012)
  *PNAS* **109**:8884–8889.
- **Siuti, Yazbek & Lu (2013)** "Synthetic circuits integrating logic and memory
  in living cells," *Nat. Biotechnol.* **31**:448–452, doi:10.1038/nbt.2510 —
  *(DOI/metadata verified via CrossRef; full text and abstract NOT retrievable this
  session — closed access).* Recombinase memory is widely reported to persist over
  many cell generations (a "~90+ generations" figure is commonly cited); that
  specific number could not be confirmed from the source in-session and is used here
  only as a literature-reported order-of-magnitude reference, NOT a session-verified
  measurement or a quantitative validation datapoint.
- **Gonzalez-Colell & Macia (2025)** "Characterization of recombinase activity
  across cellular growth phases," *Sci. Rep.* **15**:31311,
  doi:10.1038/s41598-025-17024-y — *(full text)* increased Bxb1 expression
  imposed **no measurable metabolic burden** on the control construct.
- **Bienick, Young, Klesmith, Detwiler et al. (2014)** "The Interrelationship
  between Promoter Strength, Gene Expression, and Growth Rate," *PLoS ONE*
  **9**:e109105, doi:10.1371/journal.pone.0109105 — *(full text)* quantified the
  approximately linear expression–growth tradeoff (bounds the reporter fitness cost).
- Flip kinetics inherited from Module 3 (Hsiao et al. 2016, *Mol. Syst. Biol.*
  **12**:869).
