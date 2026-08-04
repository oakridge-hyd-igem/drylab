# Section 6 Numerical Audit Appendix

_BioChronos Dry Lab. As-of 2026-08-02. Prepared in response to Dry-Lab Revision Instructions M2 v3b (Section 6)._

Every headline result below is traced to its equation, parameter set, source-data row, and exact unrounded value, with a sanity check. Each conclusion is classified **supported**, **inconclusive**, or **contradicted**.

---

## A. Time to 1% false positive (5.4, 15.0, 15.6 min)

**Code reference.** `module2_leak_model.py`: `time_to_fp(F_thresh, beta_leak, gamma)` -> `cumulative_fp()` -> `hazard(I_leak(t))`.
**Formulae.** I_leak(t) = (beta_leak/gamma)(1 - e^{-gamma t}); hazard lambda(I) = k_flip * I/(I+K_I); F(T) = 1 - exp(-integral lambda dt).
**Parameters.** k_flip=0.4/h, K_I=10 copies; per design (beta_leak, gamma): untagged (100, 0.30), LVA tag (100, 1.34), weak RBS (10, 0.30), combined (10, 1.34).

**Grid check (the audit's central question).** The published routine used a 0.01 h (0.6 min) grid with `np.argmax(F>=thresh)` and no interpolation, so every reported value was quantised to the nearest 0.6 min. Recomputing on a 1e-4 h grid with linear event interpolation:

| design | beta_leak | I_ss | coarse (published) | fine + interp | shift |
|---|---|---|---|---|---|
| untagged | 100 | 333.3 | 5.40 min | 5.329 min | -0.07 |
| LVA tag | 100 | 74.6 | 5.40 min | 5.388 min | -0.01 |
| weak RBS | 10 | 33.3 | 15.00 min | 14.631 min | -0.37 |
| combined | 10 | 7.46 | 15.60 min | 15.207 min | -0.39 |

The shifts are all < 0.4 min. **The published values are not a grid artifact**; the coarse grid rounds up by at most one grid step. (Best practice going forward: report the interpolated values.)

**Why 333 and 75 leak units give the same 5.4 min.** untagged (I_ss=333) and LVA tag (I_ss=75) share **beta_leak=100** and differ only in gamma. At early times gamma*t << 1, so I(t) ~ beta_leak * t regardless of gamma: the two accumulate integrase identically over the first few minutes, and 1% FP is reached before the steady states diverge. The controlling variable is the **production rate beta_leak (set by RBS), not the steady state I_ss (set by beta_leak/gamma).** Early-time closed form: t_1% = sqrt(2 * F * K_I / (k_flip * beta_leak)) = 4.24 min for both (the simulated 5.3-5.4 includes the saturating denominator). LVA changes I_ss 4.5-fold but leaves the early trajectory, and therefore the time, essentially unchanged.

**Why a 44-fold leak reduction gives only a 2.9-fold time extension.** untagged vs combined differ 44.7-fold in I_ss (333 -> 7.46) but only **10-fold in beta_leak** (100 -> 10); gamma accounts for the rest of the I_ss ratio. Since t_1% ~ 1/sqrt(beta_leak), a 10-fold beta_leak drop predicts sqrt(10) = 3.16-fold more time; observed fine-grid ratio 15.207/5.329 = 2.85-fold. **The 44-fold figure is the wrong comparator** for the trust window: I_ss governs the eventual plateau height, beta_leak governs how fast 1% is first crossed. Quoting the 44-fold leak reduction next to the 2.9-fold time extension conflates two different quantities.

**Classification: SUPPORTED** (values reproduce within 0.4 min; both apparent paradoxes are explained mechanistically). Recommended wording fix: report interpolated minutes and stop pairing "44x leak reduction" with the time-extension figure.

---

## B. Trust-horizon duplication (M3 untagged 0.26 h vs M2 combined 15.6 min)

**The two numbers have different definitions, designs, and thresholds:**

| quantity | source | design | threshold | code |
|---|---|---|---|---|
| M3 trust horizon | module3_results.csv row `untagged` | untagged (bl=100, g=0.30) | **5%** flip (F_detect) | `horizon()` in M3 regrounding code |
| M2 time-to-1%-FP | module2_results.csv row `combined` | combined (bl=10, g=1.34) | **1%** FP | `time_to_fp()` |

M3 untagged crosses 5% at 0.2564 h = **15.39 min**; M2 combined crosses 1% at 0.2535 h = **15.21 min**. Both round to ~15.6 / 0.26 on their respective coarse grids. They are numerically close but arise from different designs at different thresholds.

**Derived or coincidental? A perturbation test settles it.** If the equality were a shared derived quantity, both would move identically under any shared parameter. Bumping k_flip 0.4 -> 0.5: M3 untagged horizon -> 13.18 min, M2 combined t_1% -> 13.43 min. **The two cross over** (M3 was higher at baseline, M2 is higher after the bump). A derived identity cannot change sign. **The equality is coincidental**, not a copy or a shared bug: a design with 10x lower leak but a 1% (stricter) threshold happens to land near a design with 10x higher leak read at a 5% (looser) threshold.

**Classification: SUPPORTED as independent** (equality is coincidental; both trace cleanly to distinct rows and formulae). No shared-variable error. Recommended: annotate both tables to state the threshold and design so the coincidence is not read as a duplication.

---

## C. Memory-stability logic (~100% with leak vs ~60% leak-free)

**Code reference.** `module4_memory_model.py`: deterministic `dfdt(f,t,lam,c) = (1-f)(lam - MU*c*f)`; stochastic `wright_fisher(N,f0,c,lam_h,n_gen,reps)`; `f_star(lam,c)=min(1, lam/(MU*c))`. MU=0.30/h, gen time TAU=2.31 h, f0=0.90, c=2%, 90 gen = 207.9 h.

**The audit's valid objection:** "leak reinforcing an already-written record does not prove stability is non-binding after leak suppression." Correct. The original claim compared leak-ON (100%) against the artificial leak-free case (60%), which is not the operating point of a **leak-suppressed** design. Showing the leak-suppressed scenario explicitly:

| scenario | lambda (/h) | 90-gen retention | mechanism |
|---|---|---|---|
| leak ON (untagged) | 0.388 | 100% | record continuously replenished |
| **leak-SUPPRESSED (combined)** | **0.171** | **100%** (f* pins at 1.0) | residual leak still exceeds selection drain MU*c |
| TRUE leak-free (lambda=0), deterministic | 0 | 72.1% | pure selection erosion, mean-field |
| TRUE leak-free (lambda=0), Wright-Fisher N=1000 | 0 | 61.8% (published 60.4%) | selection + genetic drift, small-N worst case |

The combined design is leak-**suppressed but not leak-free**: its steady-state integrase I_ss is 44-fold lower than untagged (333 -> 7.46), but because the flip hazard saturates (lambda = k_flip * I/(I+K_I)), the **hazard lambda itself is only ~2.3-fold lower** (0.388 -> 0.171/h), not 44-fold. This is the same I_ss-vs-hazard distinction drawn in Section A: I_ss sets plateau height, not rate. Even at this modestly reduced hazard, lambda=0.171/h still exceeds the selection drain (MU*c = 0.30*0.02 = 0.006/h) by ~28-fold, so f* = min(1, lambda/(MU*c)) pins at 100%. Retention only erodes in the hypothetical **zero-leak** corner, which no real construct occupies.

**Two mechanisms the audit asked to separate:**
1. **DNA-state stability** (does the inverted cassette revert molecularly?): Bxb1 inversion is irreversible in the absence of RDF; the model has **no molecular reversion term**. The DNA record itself is permanent.
2. **Population selection** (do flipped cells get outcompeted?): the MU*c*f term. This, not DNA reversion, is what erodes the leak-free line to 60-72%. It is a population-fraction effect, not a memory-integrity effect.

**Classification: INCONCLUSIVE for the leak-free corner, SUPPORTED for realistic operation.** Memory stability is non-binding **for any design that retains even suppressed leak** (f* = 100%). It becomes potentially binding only in a true zero-leak construct under fitness cost, where population selection (not DNA reversion) erodes the flipped fraction to ~60-72% over 90 generations. The deprioritisation of memory stability holds for the current constructs but should be stated as conditional on residual leak, and DNA-state permanence should be reported separately from population-fraction retention.

---

## Conclusion table

| # | headline result | classification | one-line basis |
|---|---|---|---|
| A | 5.4 / 15.0 / 15.6 min to 1% FP | **supported** | reproduces within 0.4 min on a fine grid; beta_leak (not I_ss) sets the time |
| A | "44x leak reduction -> 2.9x time" pairing | **contradicted (as stated)** | time scales with beta_leak (10x -> sqrt10 = 3.16x), not I_ss (44x); comparator is wrong |
| B | M3 0.26 h == M2 15.6 min | **supported as coincidental** | different design/threshold; the two cross over under a k_flip perturbation |
| C | 90-gen retention ~100% (leak) vs ~60% (leak-free) | **inconclusive -> supported (conditional)** | leak-suppressed design still pins f*=100%; erosion only in the zero-leak corner, and it is population selection, not DNA reversion |

**Net:** two headline numbers are numerically sound (A, B); one framing is wrong and should be dropped (the 44x/2.9x pairing); one conclusion (C) needs the conditional "for any design retaining residual leak" qualifier and a DNA-vs-population distinction. None of the three is contradicted at the level of the underlying model output.
