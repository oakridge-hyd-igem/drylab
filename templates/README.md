# BioChronos Dry Lab: Analysis Templates for Incoming Wet-Lab Data

Prepared 2 September 2026, ahead of the 4 September M1/M2 construct delivery, per the
Rescue HQ timeline decision gates (9-14 Sep pilot data, 15-25 Sep biological replicates).

## Purpose

Each template pairs a CSV schema with a ready-to-run analysis script, so that dropping a
real plate-reader export or PCR log into the matching CSV produces fitted parameters and
updated figures the same day, without writing analysis code after data arrives.

All three scripts were self-tested against synthetic data generated from this project's own
grounded model parameters (Module 1's Hill fit, Dry Model D2's pBAD leak predictions, Module
4's leak-on memory retention) before being handed off, to confirm they recover the expected
values and produce readable output.

## Templates and scripts

| # | Template CSV | Script | Answers | Feeds into |
|---|---|---|---|---|
| 1 | M1_dose_response_template.csv | fit_M1_dose_response.py | Real Hill K and n for the copper sensor | Module 1 parameter update |
| 2 | M2_pilot_template.csv | fit_M2_pilot.py | Real OFF-state leak and induced-switching level, Untagged vs LVA | Dry Model D2 leak comparison, the LVA-vs-Untagged call |
| 3 | M2_replicates_memory_template.csv | fit_M2_replicates_memory.py | Population-to-population CV (task T3) and memory retention after inducer removal | Module 5 single- vs two-integrase decision, Module 4 memory validation |

## Usage

1. Copy the matching template CSV, delete the example rows, fill in real data.
2. Run: `python <script>.py <your_filled_csv>.csv`
3. Each script writes its own summary CSV and comparison PNG (measured vs the literature-proxy
   baseline already used in the models) into the working directory.

## Reporting-language rules baked into template 2 and 3

Per the Rescue HQ Section 6 rules: `pcr_call` must be one of `OFF_only`, `mixed`, `ON_only`.
This is an ON-sequence detection frequency category, not a per-cell or per-molecule flip
fraction. The scripts enforce this category set and will raise an error on any other value.

## Column definitions

- `cu_uM`: copper concentration in micromolar.
- `GFP_OD_normalized`: GFP_raw / OD600, background-subtracted if applicable; state your
  subtraction method in a comment row if used.
- `tag_variant`: `untagged` or `LVA` (matches M2 v3b-Untagged / M2 v3b-LVA naming).
- `induction`: `glucose`, `no_inducer`, `arabinose`, or `arabinose_then_removed` (memory test).
- `bio_replicate_n`: independent culture/colony, not a technical (pipetting) replicate.
- `hours_since_inducer_removal`: leave blank (NaN) for rows that are not part of the memory test.
- `pcr_call`: orientation-specific endpoint colony PCR result category (see above).

## Known limitation to fix once real data exists

The example rows in each template use illustrative numbers, not measurements. Delete them
before adding real data; the scripts skip lines starting with `#` but will otherwise try to
fit the example rows as if they were real.
