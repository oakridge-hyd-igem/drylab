# Dry Lab Modelling Plan
## Guiding Principle
 
Each model maps to a specific design decision and answers it through simulation, so construct and parameter choices rest on quantitative evidence rather than guesswork.
 
The system is a biological recorder: engineered *E. coli* that converts a transient contamination event into a permanent, irreversible DNA flip (via Bxb1 integrase), read out later. Because the flip is irreversible, the modelling must guard against two failure modes:
 
- **False positives** — recording an event that never happened, caused by integrase leakiness.
- **False negatives** — missing an event that did happen because it was too brief.
The plan predicts and balances these two failure modes and confirms that a recorded state persists.
 
---
 
## The Failure Modes
 
| Failure mode | Question | Model |
| --- | --- | --- |
| Records something that didn't happen | How do we stop false recordings? | Model 1 — Leak / false positive |
| Doesn't know the concentration boundary | How much contamination triggers a record? | Model 2 — Sensor threshold |
| Misses something that did happen | How short an event can still be recorded? | Model 3 — Transient capture |
| Forgets something it recorded | Does the record survive over time? | Model 4 — Irreversibility / memory stability |
| Can't tell apart levels of exposure | Are exposure levels distinguishable? | Model 5 — Multi-state separability |
 
Models 1 and 3 are in tension: stronger leak suppression lowers false positives but weakens the response to brief signals, raising false negatives. The target construct minimises false positives while still capturing short pulses.
 
---
 
## Module 1 — Sensor Input Model (deterministic ODEs)
 
**Models:** The sensor layer. Analyte binds the regulator, derepresses the promoter, and transcription and translation produce integrase. A Hill function describes the promoter's response to analyte concentration; mass-action terms describe transcription, translation, and protein decay.
 
**Delivers:** A dose-response curve, the activation threshold, and the response sharpness (Hill coefficient). Defines the meaningful concentration range and the State 0 → recording boundary.
 
**Tooling:** Python with SciPy (`solve_ivp`), or Tellurium / COPASI.
 
---
 
## Module 2 — Leakiness / False-Positive Model
 
**Maps to:** The integrase-leakiness suppression strategies (weak RBS, LVA degradation tag, combined) and the Gold Medal leak-reduced cassette part.
 
**Significance:** Leak is the main threat to an irreversible recorder. A small amount of integrase produced with no signal present will eventually flip the DNA and permanently corrupt the record. With three competing suppression strategies, the model ranks them before construction.
 
**Models:** Integrase concentration over time under no-signal conditions, for each strategy:
 
- **Baseline:** leaky transcription, integrase accumulates.
- **Weak RBS:** lower translation rate, slower accumulation.
- **Degradation tag (LVA):** higher decay rate, lower steady-state level.
- **Combined:** both effects.
The model links integrase level to flip probability. Because the flip is irreversible, the relevant quantity is the cumulative probability of flipping over the full deployment window, not the instantaneous integrase concentration. A small leak over weeks can still corrupt the record. The model predicts flip probability as a function of leak rate × deployment time and identifies which strategy keeps false-positive flips below a set threshold (e.g. <1%) over a multi-week deployment.
 
**Delivers:** A ranked recommendation of suppression strategy, and the characterisation data for the leak-reduced cassette part (e.g. predicted fold-reduction in false-positive rate of the combined design over baseline).
 
---
 
## Module 3 — Stochastic Flip Model
 
**Maps to:** The single-bit memory flip and transient exposure capture.
 
**Why stochastic:** Flipping is a discrete, probabilistic event in each cell, and the readout depends on the fraction of the population that flipped. Deterministic ODEs do not capture this. The model uses the Gillespie algorithm (stochastic simulation algorithm) or an agent-based model in which each cell flips probabilistically based on its integrase level.
 
**Answers:**
 
- **Transient capture (core claim):** Whether a short pulse (10 min, 30 min, 1 hr) produces enough integrase to flip a detectable fraction of cells. The model simulates pulses of varying duration and intensity and predicts the minimum recordable pulse. If a 10-minute pulse does not flip a detectable fraction, the sensitivity is redesigned before the construct is finalised.
- **Population threshold:** The conditions that push the flipped fraction above a usable detection threshold.
**Detection threshold (definition of "detectable"):** "Detectable fraction" is defined by the readout method, not chosen arbitrarily, and is set as a configurable parameter in the model:
 
- **Flow cytometry:** ~1–5% GFP-positive cells (resolves a small flipped subpopulation against the negative majority). Default: 5%.
- **qPCR:** more sensitive than gel PCR; threshold set against the no-flip control.
- **PCR orientation (gel) assay:** ~10–20% flipped, below which a faint flipped band is unreliable on a gel.
The exact floor depends on the negative-control background (autofluorescence spread for cytometry, baseline band for PCR) and should be confirmed against actual control data rather than assumed. The model default is 5% flipped (flow cytometry); the threshold parameter is changed to match whichever readout the wet lab uses.
 
**Priority:** Highest. It underpins the core claim and feeds Modules 4 and 5.
 
---
 
## Module 4 — Irreversibility / Memory Stability Model
 
**Maps to:** The memory-persistence claim — a written state must survive after the signal disappears.
 
**Why it matters:** A flip is only useful if it stays readable through the rest of the deployment, across cell division, dilution, and growth. The model tests whether the flipped state remains stable on that timescale.
 
**Models:** The flipped subpopulation over time under growth, testing whether the flipped fraction is preserved or eroded as the culture divides over days to weeks. The recombined DNA is genetically permanent, so the question is population-level readability over the deployment window, not reversal in individual cells.
 
**Delivers:** A stability curve over deployment time, and any timescale at which dilution or overgrowth threatens readability. Supports the "persistent, retrievable record" claim.
 
---
 
## Module 5 — Multi-State Separability Model
 
**Maps to:** The 0 / 1 / 2 multi-state system.
 
**Models:** Using the Module 3 engine, the population distribution of flipped fraction under low/transient versus high/sustained exposure, then tests whether the State 1 and State 2 distributions are separable or overlapping.
 
**Delivers:** A prediction of whether threshold tuning produces separable states or whether a two-integrase architecture is required.
 
**State 1 vs State 2 definition (open design question):** A single cell either flips or it does not, so "partial recombination" at State 1 means the fraction of the population that flipped. The cleanest definition is a population-level threshold (e.g. 20–60% flipped = State 1, >60% = State 2). The alternative is a two-stage system using two integrases with two thresholds. The two require different constructs; the model quantifies which approach yields separable states.
 
**Metric:** Distribution overlap coefficient, or classification accuracy between states.
 
---
 
## Module 6 — Deployment / Encapsulation Model
 
**Models:** Analyte diffusion into the hydrogel or filter-paper encapsulation matrix. A spike in the water reaches the bacteria with delay, and the matrix lengthens and smooths the effective pulse. A 1D reaction-diffusion model captures this.
 
**Value:** Supports the transient-capture claim, since the matrix stretches a sharp spike into a longer effective exposure, and connects the modelling to deployment conditions. Lower priority; depends on the core models being complete.
 
---
 
## Module 7 — Sensitivity Analysis
 
**Maps to:** Every model. The conclusions depend on parameter values (Bxb1 kinetics, leak rate, RBS strength, degradation rate, Hill coefficient, binding affinity) drawn from literature or estimated. This module tests how far those conclusions hold when the parameters shift.
 
**Method:** One-at-a-time (OAT) sweep. Each key parameter is varied across a plausible range (e.g. ±1 order of magnitude, or the spread reported in the literature) while the others are held fixed, and the effect on the model output is recorded.
 
**Parameters per module:**
 
- Module 1: Hill coefficient, regulator binding affinity.
- Module 2: leak rate, degradation rate.
- Module 3: integrase accumulation rate, flip threshold.
**Delivers:**
 
- **Robustness:** which conclusions hold across the whole range (the strong claims).
- **Sensitivity:** which parameters the output depends on most. These are the values worth measuring carefully in the wet lab.
**Limitation:** OAT does not capture interactions between parameters shifting together. A global method (Sobol, Latin-hypercube) would, but is out of scope here. OAT is sufficient for this project; the limitation is stated rather than hidden.
 
---
 
## Prioritisation of Modules
 
1. **Module 3 (stochastic flip)** — underpins the core claim and feeds Modules 4 and 5. Build first.
2. **Module 2 (leakiness)** — ties to the Gold Medal part and the suppression-strategy decision.
3. **Module 1 (sensor)** — defines the activation threshold.
4. **Module 4 (irreversibility)** and **Module 5 (multi-state)** — build on Module 3.
5. **Module 7 (sensitivity analysis)** — run on each model as it is completed.
6. **Module 6 (deployment)** — if time allows.
---
 
## Model Summary Cards
 
### Card 1 — Leak / False-Positive Model
 
**Question:** How do we stop the system from recording events that never happened?
 
**Inputs:**
 
- Promoter leakiness
- RBS strengths
- Protein degradation rate (with / without LVA tag)
- Literature values for Bxb1 activity
- Deployment duration
**Outputs:**
 
- Cumulative false-positive flip rate over deployment time
- Comparison of baseline, weak RBS, LVA tag, and combined designs
**Use:**
 
- Selects the leak-suppression strategy
- Reduces false recordings
- Provides quantitative support for engineering decisions
### Card 2 — Sensor Threshold Model
 
**Question:** What contamination concentration is enough to trigger a recording?
 
**Inputs:**
 
- Analyte concentration range
- Regulator binding affinity
- Promoter / Hill parameters
**Outputs:**
 
- Dose-response curve
- Activation threshold (State 0 → recording boundary)
- Response sharpness (Hill coefficient)
**Use:**
 
- Defines the meaningful concentration range
- Sets the boundary between no event and a recording
- Provides the concentration axis of the operating range (paired with Card 3's duration axis)
### Card 3 — Transient Capture Model
 
**Question:** How short can an exposure be and still be recorded?
 
**Inputs:**
 
- Exposure duration (10 min, 30 min, 1 hr, etc.)
- Integrase accumulation kinetics (build-up and decay rates)
- Literature values for gene expression
**Outputs:**
 
- Percentage of cells that flip after different exposure times
- Minimum exposure needed for reliable detection
- Operating range of the system
- Effect of leak suppression (Card 1) on the minimum recordable pulse
**Use:**
 
- Predicts whether short contamination events can be recorded
- Identifies the minimum recordable event
- Supports the core goal of recording past exposures
### Card 4 — Irreversibility / Memory Stability Model
 
**Question:** Does the recorded state survive over the deployment, or does it fade?
 
**Inputs:**
 
- Cell growth / division rate
- Flip stability assumptions (genetic permanence of the recombined state)
- Deployment duration
- Initial flipped fraction from Module 3
**Outputs:**
 
- Stability of the flipped fraction over time
- Timescale over which the record remains readable
- Any point at which dilution or overgrowth threatens readability
**Use:**
 
- Supports the persistent-record claim
- Defines how long after an event a cartridge can still be read
- Tests the irreversibility component of the system
### Card 5 — Multi-State Separability Model
 
**Question:** Can the system distinguish different levels of exposure (State 0 / 1 / 2)?
 
**Inputs:**
 
- Population flip distributions under low/transient vs high/sustained exposure (from Module 3)
- State threshold definitions (% flipped boundaries)
- Exposure intensity and duration ranges
**Outputs:**
 
- State distribution under each exposure regime
- Whether State 1 and State 2 are separable or overlapping
- Recommendation: threshold tuning vs two-integrase architecture
**Use:**
 
- Predicts whether the 0/1/2 scheme will work before construction
- Prevents committing to a multi-state design that cannot separate cleanly
- Informs the State 1 vs State 2 mechanism decision
### Card 6 — Sensitivity Analysis
 
**Question:** Do the conclusions hold when the parameters shift?
 
**Inputs:**
 
- Key parameters from each model (leak rate, degradation rate, integrase accumulation rate, flip threshold, Hill coefficient, binding affinity)
- Plausible range for each (literature spread, or ±1 order of magnitude)
**Outputs:**
 
- Which conclusions are robust across the parameter range
- Which parameters the output is most sensitive to
- Ranked list of parameters worth measuring carefully in the wet lab
**Use:**
 
- Answers the "what if your parameters are wrong?" question directly
- Identifies the strong, range-independent claims
- Directs wet-lab measurement effort toward the parameters that matter most
---
