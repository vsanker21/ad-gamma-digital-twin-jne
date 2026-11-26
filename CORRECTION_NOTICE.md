# Correction and supersession notice (v2.0.0)

**Date (UTC):** 2026-11-26  
**Current Zenodo version DOI:** `10.5281/zenodo.24570086` (v2.0.0)  
**Concept DOI:** `10.5281/zenodo.22746384`  
**Supersedes:** GitHub release `v1.0.0` and Zenodo record `10.5281/zenodo.22746385` (JNE resubmission archive).

## Why v1 must not be used

The v1 analysis-ready tables and scripts were insufficient for confirmatory claims:

1. Outcome provenance was not machine-traceable (rows labeled only `published_cohort_table`).
2. Iaccarino outcomes were min–max normalized, creating artificial ties at 0.0 and hiding distinct raw Aβ values.
3. Unpaired SWR “intrinsic gamma” priors were attached to ELISA animals.
4. Etter S1 (`n=7`) was stated in the plan but the deposited script analyzed all non-null rows (`n=10`).
5. Incomplete Etter mice (181/183/184) had 10 Hz ratios used as if they were 40 Hz responses in prior pipelines.
6. The “preregistration” JSON embedded finished coefficients and pass/fail flags.
7. Power calculations were attached to the wrong sample size for S1.
8. Phase 2 headline 28 Hz conflicted with equal-usable weighting (42 Hz) and noise bootstrap (0% at 28 Hz).

## What v2 changes

- Raw published values with exact source panel/table citations.
- Iaccarino primary outcome = fractional Aβ reduction vs Dark 1 hr VC mean (and raw pg/ml sensitivity).
- No per-animal endogenous gamma column on the Iaccarino ELISA table.
- Etter S1 inclusion coded as complete 10+40+80 only; incomplete mice retained but flagged.
- Analysis specification separated from results; no backdated fake preregistration.
- Claim boundary explicitly excludes personalization validation and pooled ML confirmatory inference.

## Claim boundary (read this)

v2 supports a **frequency-titration consistency** statement for Iaccarino assigned 20/40/80 Hz arms.  
v2 does **not** support a claim that personalized gamma entrainment was tested or validated.
