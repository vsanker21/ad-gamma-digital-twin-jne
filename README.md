# ad-gamma-digital-twin — corrected verified archive v2.0.0

**Status:** Corrected scientific archive after JNE-111725 rejection.  
**Current Zenodo version DOI:** https://doi.org/10.5281/zenodo.24570086 (v2.0.0)  
**Concept DOI (resolves to latest):** https://doi.org/10.5281/zenodo.22746384  
**Supersedes:** GitHub release v1.0.0 and Zenodo 10.5281/zenodo.22746385 — **do not use**.

Read `CORRECTION_NOTICE.md` before using any numbers.

## What this archive contains

| Path | Description |
|------|-------------|
| `data/etter_2019_moesm4_fig2n_raw.csv` | Raw stim/baseline ratios + inclusion flags |
| `data/iaccarino_2016_ext_data_table1_raw.csv` | Raw Aβ1-40 ELISA + dark-baseline fractional reduction |
| `data/ANALYSIS_SPECIFICATION.json` | Analysis rules (no results embedded) |
| `results/CORRECTED_ANALYSIS_RESULTS.json` | Statistics from the specification |
| `data/phase2_weighting_sensitivity.json` | Prior Phase-2 sensitivity (42 Hz equal-usable; 28 Hz GC-only) |
| `provenance/DATA_PROVENANCE.md` | Source ledger |
| `code/` | Rebuild + analysis scripts |

## Primary scientific statement supported by v2

In the Iaccarino 2016 VC ELISA cohort (n=18), greater assigned stimulation mismatch from 40 Hz is associated with lower fractional Aβ1-40 reduction versus the Dark 1 hr VC baseline (Spearman with permutation inference). This is a **frequency-titration consistency** result.

## Statements NOT supported by v2

- Personalized / individualized gamma entrainment validation
- Independent predictive validation on unseen experiments
- Robust 28 Hz stimulation optimum
- Confirmatory machine learning on pooled n=28

## Contact

Vivek Sanker — vsanker@stanford.edu  
Atman Desai — atman@stanford.edu
