# Module 3 — Validation against Hsiao 2016 (independent mechanistic benchmark)

**BioChronos** · Module 3 (Stochastic Flip) validation pass

## Purpose
Module 3's parameters are *anchored* to Hsiao et al. 2016 (k_flip=0.4 h⁻¹,
k_deg=0.3 h⁻¹, Kd=10 molecules — all verified verbatim from the fetched full text,
PMC/DOI 10.15252/msb.20156663). But anchoring parameters is not the same as
validating behaviour. This pass tests the one structural simplification the model
makes, against Hsiao's own mechanistic equation.

## Parameter check (no discrepancy)
Hsiao's **nominal** parameter set is explicitly "k_flipA = k_flipB = 0.4 h⁻¹,
k_deg = 0.3 h⁻¹ (2.3 h half-life), Kd = 10 molecules" (quoted from full text). Our
baseline k_flip=0.4 matches this. (The 0.2/0.3 h⁻¹ values elsewhere in the paper
are a later *asymmetric two-integrase* refinement for a specific fitting exercise,
not the base flip rate — so there is no discrepancy with our baseline.)

## The structural test: first-order vs tetramer hazard
- **Module 3 uses** a first-order saturating flip hazard: λ(I) = k_flip·I/(I+K_I).
- **Hsiao Eq. (2) uses** a 4th-order tetramer binding hazard — a flip requires four
  integrase monomers to occupy attB/attP and form a tetramer:

  λ(I) = k_flip · I(I−1)(I−2)(I−3) / [Kd⁴ + Kd³I + Kd²I(I−1) + KdI(I−1)(I−2) + I(I−1)(I−2)(I−3)]

I implemented Hsiao's exact equation from the paper and compared.

### Result — the two agree where it matters for genuine events, diverge at low integrase
| Integrase I | first-order λ (per h) | Hsiao tetramer λ (per h) |
|---|---|---|
| 5 | 0.133 | 0.003 |
| 10 | 0.200 | 0.049 |
| 20 | 0.267 | 0.184 |
| 33 | 0.307 | 0.269 |
| 50 | 0.333 | 0.316 |

The tetramer hazard has a **steep threshold below ~4 integrase** (mathematically it
is zero for I<4) and half-maximal near I≈18, versus our first-order half-max at
I=K_I=10. At high (induced) integrase the two converge.

### Impact on the recorder conclusions
False-positive flip at 24 h at each design's leak steady-state integrase:

| Design | I_leak | first-order (Module 2 used) | Hsiao tetramer |
|---|---|---|---|
| untagged | 33.3 | 99.9% | 99.9% |
| LVA tag | 7.5 | 98.4% | 35.1% |
| weak RBS | 3.3 | 90.9% | ~0% |
| combined | 0.75 | 48.8% | ~0% |

**Two conclusions:**
1. **Genuine-event recording is validated.** At induced (high) integrase both
   models give the same near-complete flip — the core "a real pulse gets recorded"
   claim (Module 3) does not depend on the hazard form.
2. **The first-order model is conservative (worst-case) for leak false positives.**
   Under Hsiao's mechanistic tetramer hazard, the leak-suppressed designs
   (weak RBS, combined) are *much safer* than Module 2 reported — the tetramer
   requires ≥4 integrase to flip, and leak keeps these designs below that. Module
   2's "leak-limited" headline is unchanged for the untagged design but the
   leak-suppressed designs gain a large safety margin under the more realistic
   mechanism.

## What this means for the module set
- Module 3's headline (recorder is leak-limited, not sensitivity-limited) **stands**
  — it was derived on the first-order hazard, which is the conservative choice.
- The first-order simplification is now **quantified, not just flagged**: it
  overestimates leak-driven false positives at low integrase and is therefore a
  safe (pessimistic) modelling choice. A mechanistic tetramer option should be a
  Module 7 structural-sensitivity variant.

## Additional Hsiao benchmarks noted (not yet reproduced)
Hsiao also reports a Δt₉₀ detection limit (~4 h for their experimental system) and a
maximum switchable population fraction of ~60% (attributed to non-fluorescent
sub-populations having a slight growth advantage). These are two-integrase
timing-gate metrics; reproducing them requires the full two-input architecture
(Module 5 territory), so they are logged as future validation targets, not run here.

## Source
- Hsiao, Hori, Rothemund & Murray (2016) "A population-based temporal logic gate
  for timing and recording chemical events," *Mol. Syst. Biol.* **12**:869,
  doi:10.15252/msb.20156663 — full text (PMC), Eq. (2) and nominal parameters
  read directly this session.
