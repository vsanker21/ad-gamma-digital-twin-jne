# Data provenance ledger (v2.0.0)

## Etter et al., 2019 (Nature Communications)

| Field | Value |
|------|------|
| Article | Etter et al. Nat Commun 2019 |
| DOI | 10.1038/s41467-019-13260-9 |
| Source file | 41467_2019_13260_MOESM4_ESM (fig_2n sheet/export) |
| Local source used | `G:\Alzheimers Dementia - Digital Twin\real_experimental_data\etter_2019\41467_2019_13260_MOESM4_ESM_fig_2n.csv` |
| n rows | 10 |
| Complete 10+40+80 | 7 |
| Incomplete (10 Hz only) | Etter_181, Etter_183, Etter_184 |
| Preferred frequency | argmax among measured stim/baseline ratios (proxy only) |
| Plaque burden | **Not deposited for confirmatory use** (anonymous % list; mouse-ID pairing unverified) |

## Iaccarino et al., 2016 (Nature)

| Field | Value |
|------|------|
| Article | Iaccarino et al. Nature 2016 |
| DOI | 10.1038/nature20587 |
| Source | NIHMS886326 Extended Data Table 1 (Aβ1-40 ELISA, VC 1 hr arms) |
| Local source used | `G:\Alzheimers Dementia - Digital Twin\real_experimental_data\iaccarino_2016\iaccarino_vc_cohort_with_endogenous_gamma.csv` |
| n rows | 18 |
| Dark 1 hr VC baseline mean (pg/ml) | 354.408332 |
| Animals with no reduction vs dark | 8 |
| Unique raw Aβ values | 17 |
| Endogenous gamma | **Excluded from animal table** (not published paired to ELISA mice) |

## Reproducibility

```bash
python code/build_verified_archive.py
python code/run_corrected_analyses.py
```

Random seed: 42
