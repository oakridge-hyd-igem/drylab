# BioChronos Dry Lab: F4 Citation-Verification Log

**As-of 2026-08-02.** Verification of the prior-art claims in the consolidated report's F4 leak-control survey, per M2 v3b Revision Instructions Section 9. Each claim is marked VERIFIED (checked against the official Registry / team page this session), CORRECTED (attribution was wrong, now fixed), or UNVERIFIED (could not be independently checked this session; requires a manual Registry/wiki check before use in a submission).

## Method and scope
- Sources checked: iGEM Registry part records via the Registry XML API (parts.igem.org/cgi/xml/part.cgi), plus web-search snippets of Registry part pages and iGEM team wikis (20xx.igem.org), retrieved 2026-08-02.
- Limitation: the Registry HTML part pages returned HTTP 403 and could not be fetched directly this session; only the XML API records and search-result snippets were available. Several XML records also returned intermittent 403 (rate-limiting). Claims that could not be reached are marked UNVERIFIED rather than assumed correct, and page-content attributions not confirmed from the actual HTML page are flagged as such.

## Claim-by-claim

### C1. "Weaker RBS on the integrase: Peking 2012 (BBa_K907000, Bxb1) noted pBAD leakage would cause flipping and placed a weaker RBS (B0033) after the promoter."
**Status: CORRECTED (attribution was wrong on two counts).**
- BBa_K907000 is a **KAIST 2012** part, not Peking. Registry XML: entered 2012-09-20, authors Dong-hui Choe and Soo-in Lee, "Mycobacterium Phage Bxb1 gp35, DNA integrase," released HQ 2013. The part page links to the KAIST_iGEM_2012 wiki.
- The 2012 KAIST characterization used a **Trc promoter with IPTG induction**, not pBAD. KAIST noted the Trc promoter had basal Bxb1 expression (color change without IPTG) and applied optimizations to reduce basal transcription. So the 2012 part did not use pBAD or the B0033 weak RBS.
- **The later pBAD/B0033/AraC characterization is attributed to Fudan China 2017 by the M2 v3b Revision Instructions (Section 9), not independently team-verified here.** What was actually retrieved from the Registry this session: the K907000 XML lists twin parts **BBa_K2243002 and BBa_K2243003** (later-season derivatives), and search-result snippets for the K907000 page mention an integrase orthogonality test using the pBAD promoter (BBa_I13453) + B0033 RBS + AraC. The team ownership of the K2243 twin parts (and of any other later-season derivatives) was **not** independently confirmed, and the K907000 HTML part page itself returned HTTP 403 and was never fetched directly. Treat the specific part-series-to-team linkage as TBD; the Peking-to-Fudan correction rests on the revision instructions' own statement plus the confirmed fact that the base part is KAIST 2012.
- **Corrected claim for F4:** do not attribute the pBAD/B0033 leakage discussion to "Peking 2012." The base part BBa_K907000 is KAIST 2012 (Trc/IPTG, documented basal integrase leakage). The pBAD/B0033/AraC integrase-leakage characterization on that part page is later work attributed to Fudan China 2017 per the revision instructions; confirm the exact team and part series on the official Registry page before final submission.

### C2. "Degradation tags: Waseda 2020 used ssrA-LVA tags."
**Status: UNVERIFIED this session (Registry/wiki 403). Requires manual check.**
- The ssrA-LVA degradation-tag concept itself is well grounded independently (Andersen et al. 1998, Appl Environ Microbiol 64:2240, ssrA-LVA t1/2 ~40 min), which is what our model uses.
- The specific attribution to "Waseda 2020" using ssrA-LVA on an integrase or reporter was NOT re-verified against the Waseda 2020 team wiki this session. Before using it in the F4 contribution claim, confirm on the official 2020.igem.org Waseda team page: which part, which protein tagged, and the reported effect.

### C3. "CFLS 2025 used a higher-efficiency LAA-LAA tag to minimise leaky expression of a toxic protein."
**Status: CORRECTED (context mismatch; not Bxb1-LVA evidence).**
- Per the revision instructions (Section 9), CFLS 2025 reports LAA-LAA degradation-tag effects in an **NhaA and toxin / suicide-switch context**, not on Bxb1 or integrase leak. It must not be presented as direct evidence for tagging Bxb1 with a degradation tag.
- **Corrected claim for F4:** cite CFLS 2025 only as general prior art that degradation tags reduce leaky expression of a burdensome/toxic protein, explicitly noting the context is NhaA / a toxin switch, not integrase memory. Do not use it to support Bxb1-LVA specifically.

### C4. "Tight multi-repressor architectures: insulated TetR/cI/LacI circuits (e.g. BBa_K2384014)."
**Status: UNVERIFIED this session (Registry 403). Requires manual check.**
- BBa_K2384014 could not be reached this session. Before using it, confirm on the official Registry page: what the part is, the team/year, and whether it is a multi-repressor insulation architecture as described.

### C5. "Choosing a lower-leak recombinase: teams comparing Bxb1/PhiC31/TP901 reported TP901 as the lowest-leak of the three."
**Status: UNVERIFIED this session. Requires a specific source.**
- This is a general field statement with no specific Registry part or team page cited. Before using it in the F4 claim, attach a concrete source (a specific team wiki or publication making the Bxb1/PhiC31/TP901 leak comparison) or reduce it to a qualitative statement with an explicit "source TBD" flag.

### C6. "None [of the surveyed teams] quantified the cumulative probability that leak alone crosses a false-positive threshold."
**Status: NOT A VERIFIABLE CLAIM AS STATED. Requires documented search scope.**
- Per Section 9, a "no prior team did X" claim cannot stand without documenting the search scope, databases, query terms, and inclusion criteria. The current survey did not do a systematic Registry search.
- **Corrected framing for F4:** state the claim as "in the parts and team wikis reviewed here (listed below), leak was treated qualitatively; we did not find a cumulative-false-positive / trust-horizon treatment," and list exactly which parts/pages were reviewed. Do not assert an absolute "no team ever" without a documented systematic search.

## Summary table
| claim | subject | status | action |
|---|---|---|---|
| C1 | BBa_K907000 weak-RBS / pBAD leak | CORRECTED (partial) | Base part = KAIST 2012 (Trc/IPTG), VERIFIED. pBAD/B0033 characterization attributed to Fudan China 2017 per revision instructions, team/part-series linkage TBD. Not Peking 2012 |
| C2 | Waseda 2020 ssrA-LVA | UNVERIFIED | check Waseda 2020 wiki; tag concept itself grounded via Andersen 1998 |
| C3 | CFLS 2025 LAA-LAA | CORRECTED | NhaA/toxin context, not Bxb1-LVA evidence; cite as general prior art only |
| C4 | BBa_K2384014 multi-repressor | UNVERIFIED | check Registry page for team/year/architecture |
| C5 | TP901 lowest-leak of three | UNVERIFIED | attach a specific source or mark qualitative/TBD |
| C6 | "none quantified cumulative leak" | NOT VERIFIABLE AS STATED | reframe as "in the pages reviewed here"; document search scope |

## Verified sources (this session)
- BBa_K907000: Registry XML (part_entered 2012-09-20; authors Dong-hui Choe, Soo-in Lee; KAIST_iGEM_2012; twin parts BBa_K2243002/003). VERIFIED from XML.
- Fudan China 2017 attribution of the pBAD/B0033/AraC characterization: stated by the M2 v3b Revision Instructions (Section 9); NOT independently verified here (K907000 HTML page returned 403; K2243/K2460 part-series team ownership not confirmed).
- Andersen et al. 1998, Appl Environ Microbiol 64:2240 (ssrA-LVA half-life): grounding for the degradation-tag model parameter. VERIFIED previously.
