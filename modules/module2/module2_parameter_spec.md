# Module 2 — Leakiness / False-Positive Model

**BioChronos** · Dry lab / modelling · Oakridge-HYD iGEM 2026

**Maps to:** the integrase-leakiness suppression strategies (weak RBS, LVA
degradation tag, combined) and the Gold Medal leak-reduced cassette part.

**Core question:** with **no signal present**, what is the cumulative probability
that leak alone flips the irreversible switch over a multi-day/week deployment,
and which suppression strategy keeps that false-positive (FP) rate below a target
(default **<1%**)?

---

## 1. Why this is a separate model from Module 3

Module 3 asked a **short-timescale, per-event** question: what is the minimum
pulse that flips ≥5% of cells within a 24 h readout. Module 2 asks the
**deployment-timescale, no-signal** question. Because the flip is **irreversible**,
the relevant quantity is the **cumulative** flip probability integrated over the
*entire* deployment window — not the instantaneous integrase level. A small leak
sustained over weeks still corrupts the record permanently.

## 2. Model

Under no signal (`A = 0`), the Hill/`β_max` production terms vanish and only leak
production and clearance remain:

```
dI/dt = β_leak − γ·I,        I(0) = 0
   ⟹  I(t) = (β_leak/γ)·(1 − e^(−γt))          # deterministic leak accumulation

λ(I) = k_flip · I / (I + K_I)                   # per-cell flip hazard (Hsiao-anchored)

F(T) = 1 − exp(−∫₀ᵀ λ(I(t)) dt)                 # cumulative false-positive fraction
```

The analytic mean-field `F(T)` was **validated against the Module 3 agent-based
stochastic engine** (leak-only, `t_pulse=0`): agreement ≤0.008 at 24 h across all
four designs, so the closed form is used for the multi-week sweep.

### 2.1 How the designs act on leak

| Design | Lever | β_leak (copies·cell⁻¹·h⁻¹) | γ (h⁻¹) |
|---|---|---|---|
| untagged (baseline) | — | 100 | 0.30 |
| weak RBS | lower leak translation (×0.1) | 10 | 0.30 |
| LVA tag | faster clearance (+γ_LVA) | 100 | 1.34 |
| combined | both | 10 | 1.34 |

The flip kinetics (`k_flip=0.4`, `K_I=10`, `γ_dil=0.3`, `γ_LVA=1.04`) are inherited
unchanged from **Module 3 (Hsiao et al. 2016-anchored)** so the two modules stay
numerically consistent on the flip mechanism. **`β_leak0=100` is the T1-grounded
copper-promoter baseline** (10% of β_max; measured natural CueR fold-change 7–18×,
Fu 2024) — this is the one parameter where Module 2 now uses a copper-specific
value rather than the Hsiao integrase-construct leak. `β_max` and the Hill terms do
not enter (no-signal analysis).

## 3. Results

Steady-state leak integrase `I_ss = β_leak/γ` is the false-positive pressure;
`λ(I_ss)` is the plateau flip hazard; **time-to-1%-FP** is the trustworthy
deployment window. **β_leak baseline is grounded in measured copper-promoter
data (T1): 10% of β_max, from the natural CueR/P(copA) fold-change of 7–18×
reported by Fu 2024 → leak fraction 5.6–14%.** (The earlier 1% value was
borrowed from Hsiao's integrase construct and is optimistic for a copper promoter;
it is retained only as the engineered-best-case floor of the Module 7 sweep.)

| Design | I_ss (copies) | I_ss reduction | λ_ss (h⁻¹) | rate reduction | **time to 1% FP** | FP @ 24 h |
|---|---|---|---|---|---|---|
| untagged | 333.3 | 1× | 0.388 | 1× | **5.4 min** | 100.0% |
| LVA tag | 74.6 | 4.5× | 0.353 | 1.1× | **5.4 min** | 100.0% |
| weak RBS | 33.3 | 10× | 0.308 | 1.3× | **15 min** | 99.9% |
| **combined** | **7.46** | **45×** | **0.171** | **2.3×** | **15.6 min** | **98.2%** |

### 3.1 Headline findings

1. **The recorder is leak-limited — and grounding makes it worse.** At the
   grounded copper-promoter leak rate (10% of β_max), *every* design crosses the
   1% FP threshold within ~15 minutes. Grounding β_leak in real data (10× higher
   than the old 1% assumption) **shortened every trust window ~3×** and confirmed
   leak-driven false positives are the binding constraint, not a modelling
   artifact of an optimistic parameter.

2. **Suppression still helps and ranks cleanly:**
   **combined > weak RBS > LVA tag > untagged.** The combined design cuts
   steady-state leak integrase **45×** — but at the grounded leak level this only
   extends the trustworthy window from 5.4 → 15.6 min.

3. **Weak RBS beats the LVA tag for leak.** Lowering *production* attacks
   steady-state leak integrase directly (`I_ss ∝ β_leak`), whereas faster
   *clearance* lowers `I_ss` less effectively at these levels. The two compound
   multiplicatively when combined.

4. **Suppression alone is nowhere near sufficient for multi-day deployment.** To
   hold FP < 1% the steady-state leak integrase must be **< 1.5×10⁻³ copies/cell
   over a week** (< 3.7×10⁻⁴ over 4 weeks). At the grounded baseline the combined
   design (7.46 copies) sits **~5,000× above the 1-week requirement** and
   **~20,000× above the 4-week requirement** — 10× worse than under the old 1%
   assumption. Closing that gap needs a **promoter-level** lever, not just RBS +
   tag: a tighter/repressed promoter (much lower `β_leak`), CueR overexpression
   (Pang 2020: overexpressing CueR reduced bioreporter background signal, at some
   added metabolic burden), or a copA/cus knockout — plus possibly
   a higher flip threshold `K_I` or a lower-`k_flip` construct.

### 3.2 Gold Medal part characterisation

The leak-reduced cassette (combined weak-RBS + LVA design) delivers, relative to
the untagged baseline: **45× lower steady-state leak integrase**, **~2.3× lower
false-positive flip rate at the grounded leak level**, and **~3× longer
trustworthy window** — the quantitative fold-reduction figures for the part's
characterisation. Note the between-design fold-reductions are set by the
suppression strategies and are independent of the absolute β_leak; grounding
β_leak shifts every design's *absolute* leak up 10× without changing their
*relative* ranking.

## 4. Caveats

1. **Baseline now grounded in copper data (T1 resolved).** `β_leak` is the
   T1-grounded copper-promoter baseline (10% of β_max; natural CueR/P(copA)
   fold-change 7–18×, Fu 2024). The flip kinetics `k_flip`/`K_I` remain
   Hsiao-anchored and are the next least-certain parameters (Module 3 verification
   memo). **Module 7 sweeps `β_leak` over [10, 500] copies·cell⁻¹·h⁻¹ (1%
   engineered-best-case floor → 50% worst measured natural leak)** — the
   deployment-window conclusion scales with it, and the qualitative result
   (leak-limited) is robust across this entire range (it only gets worse as leak
   rises). The single authoritative value will be wet-lab EXP1 (uninduced flip
   rate, plan-decision C1); see T1_resolution_beta_leak.md.

2. **Deterministic leak trajectory / no extrinsic noise.** `I(t)` is a single
   shared trajectory; cell-to-cell variation in leak *expression* would broaden
   the FP distribution. The mean-field `F(T)` matches the agent-based engine at
   the population level (validated), but the spread is not modelled here.

3. **Weak-RBS factor is a design assumption.** `f_rbs = 0.1` (10× weaker relative
   to B0034, the medium-strong integrase RBS in decision D4). The RBS choice
   scales basal and induced expression equally, so it changes absolute `β_leak`
   but NOT the leak *fraction* set by the promoter. True fold-change depends on the
   RBS chosen (Anderson/Salis RBS Calculator); swept in Module 7.

## 5. References

Inherited from Module 3 (CrossRef-verified):

- Hsiao, Hori, Rothemund & Murray (2016) "A population-based temporal logic gate
  for timing and recording chemical events," *Mol. Syst. Biol.* **12**:869,
  doi:10.15252/msb.20156663 — flip kinetics: k_flip = 0.4 h⁻¹, k_deg = 0.3 h⁻¹,
  Kd = 10 molecules. (Its 1–2% leak was the *previous* β_leak baseline, now
  demoted to the engineered-best-case floor of the Module 7 sweep.)

**β_leak grounding (T1, this project):**

- Fu, Li, Wang, Wang & Fang (2024) "Development of a two component system based
  biosensor...", *Commun. Biol.* **7**:1407, doi:10.1038/s42003-024-07112-6 (full
  text) — natural CueR/P(copA) copper-promoter fold-change 7–18× ⟹ leak fraction
  5.6–14%; baseline 10% (β_leak = 100). Also Pang et al. (2020) *Front. Microbiol.*
  **10**:3031 (full text) — overexpressing CueR reduced whole-cell bioreporter
  background signal (a promoter-level mitigation lever; the paper notes CueR
  overproduction adds metabolic burden, so it is a trade-off).
- J. B. Andersen et al. (1998) "New unstable variants of green fluorescent
  protein…," *Appl. Environ. Microbiol.* **64**:2240–2246 — ssrA LVA tag
  half-life (~40 min ⟹ γ_LVA ≈ 1.04 h⁻¹).
- U. Alon, *An Introduction to Systems Biology*, Ch. 1–2 — baseline production,
  dilution (γ = ln2/τ).
- Zhang, Azarin, Sarkar et al. (2022) *Nat. Commun.* **13**:4152,
  doi:10.1038/s41467-022-31538-3 — Bxb1 rate tunable over orders of magnitude
  (basis for the wide k_flip sweep range in Module 7).
