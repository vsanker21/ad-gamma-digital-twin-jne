# Manuscript claim corrections required for scientific accuracy

Use these as the only defendable claims in any resubmission or new submission.

## Title / framing

**Avoid:** “Toward personalized gamma entrainment … pre-registered frequency-mismatch analysis” as if personalization was tested.

**Prefer:** “A population-level hippocampal network model and frequency-titration consistency analysis of published preclinical 40 Hz stimulation data”

## Abstract / Significance — keep

- Population-calibrated Brian2 / network model (simulation; hypothesis-generating).
- Iaccarino VC ELISA (n=18): assigned 20/40/80 Hz titration consistency with Aβ1-40 reduction vs dark baseline.
- Explicit limitation: personalization not tested; no paired endogenous gamma for ELISA mice.

## Abstract / Significance — remove or demote

- Any implication that personalized stimulation parameters were derived from measured intrinsic gamma and then validated.
- 28 Hz as a robust / preferred stimulation optimum (equal-usable optimum is ~42 Hz; 28 Hz is GC-only and noise-unstable).
- Pooled n=28 ML as evidence of predictive utility.
- “Independent validation” language for Phase 3 (use “consistency assessment against published titration arms”).
- Fake/backdated preregistration language. Call it an analysis specification for corrected reanalysis if discussing v2.

## Results numbers to cite from v2

### P1 (Iaccarino, n=18; fractional Aβ reduction vs dark)

- Spearman ρ = −0.768
- Permutation p ≈ 3.2×10⁻⁴
- Bootstrap 95% CI ≈ [−0.901, −0.400]
- Close (40 Hz) vs far (20/80 Hz) Mann–Whitney p ≈ 5.0×10⁻⁴

### S1 (Etter complete 10+40+80, n=7; raw 40 Hz ratio)

- Included IDs: 177, 236, 238, 255, 333, 334, 354
- Excluded (10 Hz only): 181, 183, 184
- Spearman ρ = −0.577
- Permutation p ≈ 0.227 (not significant)
- Bootstrap 95% CI ≈ [−0.907, +0.228] (includes null)
- Approx power n=7: ~20% at |ρ|=0.5; ~41% at |ρ|=0.7

### Phase 2

- Report equal-usable (GC+PAC z-scored) optimum **42 Hz**.
- Disclose GC-only **28 Hz** only as sensitivity.
- Disclose two composite terms were identically zero.
- Disclose 1% noise bootstrap fraction at 28 Hz = 0.

## Data availability statement (replacement)

All materials for the corrected analyses are in GitHub release v2.0.0  
https://github.com/vsanker21/ad-gamma-digital-twin-jne/releases/tag/v2.0.0  
and Zenodo v2.0.0: https://doi.org/10.5281/zenodo.24570086  
(concept DOI https://doi.org/10.5281/zenodo.22746384).  
Zenodo v1.0.0 (10.5281/zenodo.22746385) is superseded and must not be used.

## References

Remove or verify Ref. 13 (“He J … Alzheimers Dement. 2023…”) — reviewer could not locate it; do not cite unverifiable references.
