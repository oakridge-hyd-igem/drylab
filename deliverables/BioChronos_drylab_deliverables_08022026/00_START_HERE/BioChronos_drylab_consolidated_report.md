# BioChronos Dry Lab: Consolidated Report for Prof. Nishida and the Wet Lab

**Oakridge-HYD iGEM 2026 | prepared by the dry lab**
**Revision 2 (as-of 2026-08-02). Supersedes Revision 1 (prepared 2026-07-28).**
**Target: copper (Cu2+), confirmed by the 25 July 2026 Current Project Brief.**

This document has four parts:

1. **Answers, question by question** to Prof. Nishida's requirements doc.
2. **Results summary** from the six completed models.
3. **What we need from the wet lab**, ranked.
4. **Gaps and open items.**

## Change log (Revision 1 -> Revision 2)
This revision incorporates the M2 v3b Revision Instructions (2026-08-02). Key corrections:
- **Model labels re-mapped to the current constructs.** The primary comparison is now weak-RBS untagged (M2 v3b-Untagged) vs the same weak RBS plus LVA (M2 v3b-LVA). The report's old "combined" and "weak RBS" conditions are provisional proxies only; the old "untagged" condition was strong-RBS and is not a current construct. See the separate model-to-wet-lab mapping table.
- **Headline numbers marked as model-conditional predictions**, not measured performance of the current sequences.
- **Current M2 leak re-grounded to pBAD**, not the copper promoter. Copper-promoter leak (Fu 2024) applies only to the future copper recorder (Wet M3).
- **F4 citation attributions corrected**: BBa_K907000 is KAIST 2012 (Trc/IPTG), not Peking 2012; several prior-art claims are flagged unverified. See the separate F4 citation-verification log.
- **Numerical-audit corrections applied** (Section 6 audit): the leak-vs-time scaling, the trust-horizon coincidence, and the memory-stability logic are corrected inline below.

These eight companion deliverables carry the detail: D2_remap_comparison, D2_model_wetlab_mapping, parameter_provenance, D2_remap_predictions, section6_audit_appendix, F4_citation_verification_log, claim_status_rollup, and (deferred) the next-iteration proposal.

---

# PART 1 - Answers to Prof. Nishida's questions

## Direct questions (email, 2026-06-04)

### Q1. Achievement plan for each Gold Medal criterion
**Partial (dry lab supports; full plan is team-level).** The dry lab contributes the
evidence base for the "model" and "contribution" criteria: six quantitative models, each
making a testable wet-lab prediction (Part 2), and a characterised leak-reduced cassette
for a Registry contribution. The criterion-by-criterion team plan sits above the dry lab;
we supply the modelling evidence it cites.

### Q2. Whether the new parts qualify as novel
**Wet-lab / registration-owned; dry lab supports only.** Novelty determination is not a
modelling task: an iGEM part qualifies as new by (1) sequence/composition genuinely absent
from the Registry and (2) characterisation with real experimental data. Both need inputs
the dry lab does not hold: the **final DNA sequence** (still "confirm final sequence" per
the copper brief, so the sequence-novelty check cannot be run yet) and **wet-lab
characterisation data** (EXP1/EXP3), which modelling supports but does not substitute for.

The dry lab's contribution is limited to two supporting pieces: (a) a **prior-art scan** of
existing integrase-memory and leak-control parts in the Registry, to establish what already
exists before novelty is claimed (partly covered by the leak-control survey in F4), and (b)
the **modelling evidence** for the part's datasheet (the quantitative false-positive-vs-time
behaviour of the leak-reduced cassette, which prior integrase-memory parts were not
characterised for). The novelty call itself belongs to the wet lab and the parts-registration
lead once the sequence is final and characterisation data exist.

### Q3. Anti-leak / leakiness countermeasures
**Answered, with a re-grounding note.** Dry Model D2 ranks the suppression strategies (weak
RBS, LVA degradation tag) as a quantitative false-positive model. The **current M2 comparison
is a pBAD/AraC memory test**, so its basal-leak parameter is pBAD-specific (Guzman 1995
induced/repressed ratio), NOT the copper-promoter value. The copper-promoter leak measurement
(Fu et al. 2024) is reserved for the future copper recorder (Wet M3). See feedback-item 4 for
the past-iGEM survey and how BioChronos improves, and the sensitivity analysis (Part 2) which
confirms leak is the single governing parameter of the device.

### Q4. Reconsideration of the sensor target substance
**Resolved.** The 25 July Current Project Brief locks **copper (Cu2+)**. Module 1 delivers
the copper dose-response and detection threshold on the CueR/P(copA) sensor. (The specific
CueR parts are still "confirm final sequence" per the brief, which is why the sensor K and
Hill n are flagged for the EXP3 fit.)

### Q5. Member information
**Not a dry-lab item** (team).

### Q6. A concrete experiment plan
**Supported, not owned.** The dry lab supplies the measurement priorities and design
directives (Part 3); the wet lab owns the protocol, GMO application, and safety sign-off.

## Feedback digest (Wet Lab Experiment Preparation, Notion)

### F1. Justify BioChronos through the unique value of integrase, not cumulative exposure
A normal sensor already integrates cumulative exposure, so cumulative dose alone does not
justify the project. The justification is that the integrase **permanently records a
transient event and lets you prove afterward that it happened, even after the copper is
gone.** Two models make this quantitative:
- **Module 3** sets the *minimum recordable pulse*: because the flip is irreversible and
  integrase persists ~1/gamma after the signal ends, a brief copper spike is integrated
  into a permanent DNA state and stays readable hours to days later.
- **Module 4** shows the record is a heritable DNA inversion, retained ~100% over ~90
  generations, so the event stays provable across the deployment window.

BioChronos is not "a copper sensor." It is a copper **event recorder** whose output is a
permanent, later-verifiable molecular record of a transient exposure.

### F2. A real-world scenario where transience matters
Copper contamination in water is frequently **intermittent**: corrosion of copper plumbing
releases Cu2+ pulses when water sits stagnant overnight or over a weekend and then flushes
away when flow resumes; industrial or agricultural runoff produces short copper excursions
tied to rain or discharge cycles. A grab sample or a live sensor read during normal flow
misses these, yet the transient spike is the compliance-relevant event.

**Scenario:** a distributed water-quality monitor left in a pipe or reservoir for days to
weeks. A copper spike at 2 a.m. flips the DNA; a technician sampling days later runs a PCR
orientation assay and proves a copper excursion occurred, even though current copper is
normal. This is where a steady-state sensor fails and the recorder wins.

### F3. Advantage over real-time fluorescent sensors and electrochemical loggers

| Property | Real-time fluorescent sensor | Electrochemical logger | BioChronos (integrase recorder) |
|---|---|---|---|
| Captures a transient spike when unattended | No (observed live) | Only if sampling at that instant | **Yes: permanent DNA flip** |
| Record persists after signal is gone | No (decays) | Depends on device memory/power | **Yes: irreversible, heritable** |
| Needs continuous power / electronics | Often | Yes | **No: self-contained cells** |
| Deployable as many cheap distributed units | Limited by reader cost | Limited by device cost | **Yes: culture is the sensor** |
| Read out later by a simple assay | Live optics | Download from device | **Yes: PCR / fluorescence at collection** |

**Limits:** BioChronos gives a *record that an event occurred* (and, with
multi-state, a coarse magnitude band), not a continuous concentration-vs-time trace;
readout is end-point, not live; and it has a leak-set trust horizon (Part 2, Module 2). It
complements, not replaces, a real-time sensor. Its niche is unattended, distributed,
long-deployment detection of transient events.

### F4. Survey past iGEM leak-control strategies, then show how BioChronos improves
*(from Registry parts and team wikis; attributions verified 2026-08-02, see the F4 citation-verification log for the per-claim status.)*

**What past iGEM teams did:**
1. **Weaker RBS on the integrase** - the base Bxb1 coding part BBa_K907000 is a KAIST 2012
   part (Trc/IPTG promoter, with documented basal integrase leakage). The pBAD/B0033/AraC
   integrase-leakage characterization on that part page is a later 2017 study (attributed to
   Fudan China 2017 per the revision instructions; the exact team/part-series linkage is not
   yet independently confirmed). It is NOT "Peking 2012."
2. **Degradation tags** - ssrA-LVA tags are a standard leak-control tool (Andersen 1998).
   A specific "Waseda 2020" attribution is not yet verified against the team wiki. CFLS 2025
   reports a higher-efficiency LAA-LAA tag, but in an NhaA / toxin-switch context, not on an
   integrase; cite it as general prior art only, not as Bxb1-LVA evidence.
3. **Tight multi-repressor architectures** - insulated TetR/cI/LacI circuits (e.g.
   BBa_K2384014, attribution not yet verified this session).
4. **Choosing a lower-leak recombinase** - Bxb1/PhiC31/TP901 comparisons have been reported,
   with TP901 sometimes cited as lower-leak; this needs a specific source before use.

**What they share: leak was treated qualitatively** in the parts and wikis reviewed here -
observe leak, apply a part fix, show the OFF state looks better. We did not find a cumulative-
false-positive / trust-horizon treatment in those pages (this is not an exhaustive Registry
search; see the F4 log).

**How BioChronos improves:**
- Converts leak into a **quantitative false-positive model** with a **trust horizon** (the
  deployment time before leak alone crosses a false-positive threshold) - a design metric
  the switch-based teams did not produce.
- **Grounds the leak parameter in measured data**: pBAD induced/repressed ratio (Guzman 1995)
  for the current M2 memory test, and CueR/P(copA) fold-change (Fu 2024) for the future
  copper recorder.
- **Quantifies the combined benefit** of the two most-proven fixes (weak RBS + LVA tag) and
  shows honestly that part-level suppression alone is not enough for a leak-tight week-long
  deployment - a promoter-level fix (CueR overexpression or a copA/cus-knockout host) is
  also required for the future copper recorder.
- Inherits the field's recombinase lesson: Bxb1 for unidirectionality, with TP901 as a
  documented low-leak fallback if leak proves limiting.

### F5. Respect the hard timing constraints
Acknowledged: GMO application submitted >= 2 weeks before experiments; DNA synthesis lead
time >= 2 weeks. The modelling is sequenced to deliver the design-driving answers (which
construct to build, what to measure first) before the DNA-order deadline.

---

# PART 2 - Results summary (six completed models)

All models share one flip engine (Hsiao 2016 kinetics). Numbers below are read directly from
the saved result tables.

**Read every number below as a model-conditional prediction, not measured performance of the
current M2 v3b sequences.** The model labels are the report's original ones; they map to the
current constructs only as proxies: "combined" -> proxy for M2 v3b-LVA, "weak RBS" -> proxy for
M2 v3b-Untagged, and "untagged" was strong-RBS (renamed strong-RBS untagged model, not a
current construct). The primary decision comparison is weak-RBS untagged vs weak-RBS + LVA,
delivered separately as Dry Model D2 (re-mapped, pBAD-grounded). The design directives in Part 3
are provisional pending Wet EXP-2.

### Module 1 - Copper sensor dose-response and recording threshold
Net signal-over-background detection threshold (combined design):

| Exposure | Detection threshold | Leak background |
|---|---|---|
| 1 h | 0.58 µM | 9% |
| **2 h** | **0.25 µM** | 22% |
| 3 h | 0.18 µM | 34% |

The predicted **2 h detection threshold (0.25 µM)** is the same figure as the WMC-007 copper
reporter LOD (Pang et al. 2020) and sits well below the WHO (30 µM) and China (15 µM)
drinking-water limits. Important caveat (Section 9): the WMC-007 LOD was measured with a
**chromosomal copAp::gfpmut2 reporter in a delta-copA delta-cueO delta-cusA strain on a 5 h
assay**, so it is external context, NOT validation of a 2 h threshold for the current DH5alpha
plasmid construct. This threshold is a model prediction to be tested (Wet M1). The untagged
design is unusable beyond 2 h (leak floods the readout).

### Module 2 - Leakiness / false-positive (the leak-countermeasure ranking)
Time for leak alone to reach 1% false-positive, and shortfall vs a <1% week-long deployment:

| Design | Steady-state leak integrase | Time to 1% FP | Gap vs 1-week target |
|---|---|---|---|
| untagged | 333 | 5.4 min | 222,845x |
| weak RBS | 33 | 15.0 min | 22,285x |
| LVA tag | 75 | 5.4 min | 49,891x |
| **combined** | **7.5** | **15.6 min** | **4,989x** |

Ranking: combined is best, but part-level suppression alone is ~5,000x short of a leak-tight
week-long deployment. **Leak control is the binding constraint** and needs a promoter-level
fix on top of RBS+tag. Note (Section 6 audit A): the combined design has ~44x lower
steady-state leak *level* than untagged, but only ~2.9x longer time-to-threshold, because
time-to-1%-FP scales as ~1/sqrt(leak-driven hazard), not linearly with the steady-state leak
level. The earlier "44x leak -> 44x time" reading was a units error and is corrected here.

### Module 3 - Stochastic flip (minimum recordable pulse)
Minimum recordable pulse (net over background, 0.5 h readout) and trust horizon:

| Design | Trust horizon | Min recordable pulse |
|---|---|---|
| untagged | 0.26 h | (leak floods; unusable) |
| combined | 0.67 h | ~2.2 min |

The recorder is **leak-limited, not sensitivity-limited**: because the flip integrates
signal, even sub-minute pulses are recordable if read before leak alone crosses threshold.
Validated against Hsiao's mechanistic tetramer model: our first-order hazard is the
*conservative* choice (leak-suppressed designs are even safer under the tetramer form).
Note (Section 6 audit B): the Module 3 untagged trust horizon (0.26 h) and the Module 2
combined time-to-1%-FP (15.6 min) are numerically close but are **different metrics with
different definitions**; traced independently to code, the near-equality is **coincidental**,
not a derived identity.

### Module 4 - Irreversibility / memory stability
Retention of a written record (post-event f0 = 0.90), all designs:

| Metric | Value |
|---|---|
| 48 h hold | ~100% |
| 90-generation retention (2% cost, leak on) | ~100% |
| 90-generation retention, leak-free worst case | 60% |

Two distinct effects must be separated (Section 6 audit C). The **DNA-state stability** of an
inverted record is high on its own: an irreversible flip is not reset, so the written state
persists. The **~100% (leak on) vs 60% (leak-free worst case)** gap is a *population-selection*
effect: with leak on, ongoing flipping keeps replacing any cells that lost the state, which
masks selection against burdened cells; with leak suppressed, that masking is removed and the
leak-free case is the honest worst case. So "leak reinforces the record" does not by itself
prove stability stays non-binding *after* successful leak suppression. Reading: memory
stability is a **low-risk** constraint (DNA state is stable; 60% worst-case retention over 90
generations is still usable), but it is not proven irrelevant once leak is suppressed, and the
leak-free scenario is the one to watch. Module 7 confirms the retention metric is insensitive
to every parameter swing.

### Module 5 - Multi-state separability
3-state (none / low-transient / high-sustained) classification accuracy at 1.25 h readout,
grounded expression noise (CV=0.3):

| Design | Accuracy |
|---|---|
| untagged (1 integrase) | 79% |
| **combined (1 integrase)** | **94%** |
| combined (2 integrase) | 98% |

A single leak-suppressed integrase resolves three states at ~94%; a second orthogonal
integrase adds ~4 points. **A second integrase is an optional upgrade, not a requirement**
for the binary recorder.

### Module 7 - Sensitivity analysis (which parameters matter)
OAT sweep of 10 parameters across their plausible ranges. Parameter priority:

- **Tier 1 (measure precisely, defend):** basal leak / uninduced flip rate (dominates
  trust window and threshold), Bxb1 flip kinetics (k_flip, K_I), copper sensor K and Hill n.
- **Tier 2 (confirm):** expression-noise CV (only affects the multi-state readout).
- **Tier 3 (robust, no precise value needed):** fitness cost, growth rate, LVA-tag strength,
  weak-RBS factor.

**90-generation memory retention is insensitive to every parameter (~0% swing)** - the
strongest robustness result in the project.

---

# PART 3 - What we need from the wet lab (ranked)

Ordered by how much each measurement changes the dry-lab conclusions.

### 1. Wet EXP-2 - uninduced (leak) flip level  [single most important number]
GFP/OD600 and orientation-specific endpoint colony PCR for M2 v3b-Untagged and M2 v3b-LVA in
the **uninduced (no arabinose)** condition. For M2, uninduced means no arabinose, not no
copper. This is the top-ranked driver in the sensitivity analysis; it becomes the authoritative
pBAD basal-leak value for the M2 comparison. Report as ON-sequence detection frequency and
state categories (OFF only / mixed / ON only), not a molecular OFF/ON ratio.

### 2. Wet EXP-2 - induced switching and growth  [decides the LVA call]
GFP/OD600 and OD600 growth curves for both variants under glucose repression, no inducer, and
arabinose induction. Confirms LVA lowers OFF-state leak WITHOUT unacceptably lowering induced
switching or growth. LVA is not adopted on lower leakage alone.

### 3. Wet M1 - copper dose-response curve (sensor)
GFP vs [Cu] at a concentration ladder at a fixed readout. Fits the two sensor parameters we
currently infer: half-max **K** and Hill **n**. Prediction to test: detection near **0.25 µM
at 2 h** (model-conditional; see the Wet M1 vs WMC-007 caveat in Module 1).

### 4. Population-to-population CV of the flipped fraction  [task T3]
The spread (SD or CV) of flipped fraction across >=3 replicate cultures at fixed conditions.
Confirms the multi-state verdict. Decision hinge: if measured **CV > ~0.5**, the
single-integrase 3-state readout weakens and the two-integrase architecture matters more.

### 5. Part datasheet values
B0033 weak-RBS relative strength and the pBAD basal-activity figure. Replaces the weak-RBS
factor f_rbs = 0.1 design assumption with a real number.

## Modelling-based recommendations for the wet lab to weigh
These are recommendations from the modelling, offered for the wet lab and Prof. Nishida to
consider; they are not decisions and not authorizations. Per the revision instructions
(Sections 2 and 8), no host, genotype, promoter, or RBS change is made before Wet EXP-2, and
Prof. Nishida approves any acceptance thresholds and any host change.
- **The two-construct comparison to run now is M2 v3b-Untagged vs M2 v3b-LVA** (weak RBS,
  differing only in the LVA degradation tag). The model predicts LVA lowers OFF-state leak;
  it may also lower induced switching, so LVA is not assumed superior. Wet EXP-2 decides.
- **A promoter-level leak fix** (CueR overexpression, or a copA/cus-knockout host) is a
  **future option for the copper recorder (Wet M3)**, requiring Prof. Nishida's approval and
  GMO-application review. The current DH5alpha host is unchanged for the M2 comparison.
- **Read out at ~1.25-2 h** after induction for the best transient detection and state
  discrimination (model prediction).
- **Memory stability is low-risk but not dismissed.** Per audit C, DNA-state stability is
  high, but the leak-free worst case (60% over 90 generations) is the honest scenario once
  leak is suppressed; do not deprioritize it to zero.

## Open design decision (needs Wet EXP-2 / EXP3 data)
- **What distinguishes State 1 from State 2.** A single leak-suppressed integrase already
  gives ~94% 3-state separation in the model, so a second integrase is a possible upgrade,
  not a requirement. The call depends on the measured population-to-population noise CV
  (wet-lab request 3).

---

# PART 4 - Gaps and open items

- **No wet-lab measurement of any BioChronos construct exists yet.** Every headline number is
  a model-conditional prediction; none is measured or fitted. See the parameter provenance
  table (all 18 parameters are literature-proxy or assumed).
- **Bxb1 flip kinetics are a literature proxy (Hsiao 2016), not literature-complete.** The
  true rate needs pilot wet-lab data (Wet EXP-2); the dry model is not presented as
  literature-complete, matching Prof. Nishida's own gap note.
- **Sensor K and Hill n are inferred**, not fitted to a copper dose-response - pinned by Wet M1.
- **Q1 (Gold Medal plan)** is team-level; a criterion-by-criterion evidence mapping from the
  dry lab is not yet written. **Q2 (novel-parts)** is wet-lab / registration-owned (needs the
  final sequence and characterisation data); the dry lab's role is a Registry prior-art scan
  and the part's modelling datasheet, not the novelty determination itself.
- **Module 6 (deployment / encapsulation) has not been built.** It is the lowest-priority,
  time-permitting model; the six completed models (Dry Models 1-5, 7) do not depend on it.
- **Target locked = copper** (25 July brief), but the specific CueR parts are still
  "confirm final sequence."
- Open tasks: **T3** (population-to-population noise CV, needs Wet EXP-2/EXP3). T1 (leak
  grounding) and T2 (Hill n) are resolved.
- **Next-iteration proposal is deliberately deferred** until Wet EXP-2 data exist, per the
  revision instructions (one additional design change at a time, after the current
  two-construct comparison).

---

## Sources and artifacts

**Model results**
- BioChronos Dry Models 1-5 and 7 (this project's dry-lab artifacts; Module 6 not yet built)

**Companion deliverables (M2 v3b revision, 2026-08-02)**
- D2_remap_comparison - M2 v3b-Untagged vs v3b-LVA one-page comparison
- D2_model_wetlab_mapping - corrected model-to-wet-lab mapping table
- parameter_provenance - parameter table with provenance/units/context/status
- D2_remap_predictions - prediction plots, 3 conditions x 2 variants
- section6_audit_appendix - numerical-audit appendix (Section 6 issues A/B/C)
- F4_citation_verification_log - citation-verification log for the F4 prior-art claims
- claim_status_rollup - supported / inconclusive / contradicted conclusion table
- dnafree_week_deliverables - DNA-free copper-tolerance screen (Wet-lab week)

**Grounding papers**
- Fu et al. 2024, *Commun. Biol.* 7:1407 - copper-promoter leak
- Hsiao et al. 2016, *Mol. Syst. Biol.* 12:869 - flip kinetics
- Andersen et al. 1998, *Appl. Environ. Microbiol.* 64:2240 - ssrA-LVA degradation tag
- Taniguchi et al. 2010, *Science* 329:533 - expression noise
- Pang et al. 2020; Stoyanov et al. 2001 - copper sensor
- Bienick 2014; Gonzalez-Colell & Macia 2025 - fitness cost

**iGEM leak-control precedents** (see F4 citation-verification log for per-claim status)
- Registry BBa_K907000 (KAIST 2012, Trc/IPTG; pBAD/B0033 characterization attributed to
  Fudan China 2017, team/part-series linkage TBD)
- Registry BBa_K2384014 (attribution unverified this session)
- Team wikis: Waseda 2020 (unverified), CFLS 2025 (LAA-LAA in NhaA/toxin context)
- Guzman et al. 1995, *J. Bacteriol.* 177:4121 - pBAD induced/repressed ratio (current M2 leak)

**Project brief**
- 25 July 2026 Current Project Brief (copper)
