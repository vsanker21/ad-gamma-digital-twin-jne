"""
Build scientifically corrected Phase-3 tables, analysis specification, and results.

Design principles (post JNE-111725 rejection):
1. Deposit RAW published values with exact source citations.
2. Do not attach unpaired endogenous-gamma priors to Iaccarino ELISA animals.
3. Etter S1 inclusion is coded as complete 10+40+80 Hz measurements only.
4. Analysis specification is separate from results (no results in the plan file).
5. Claims are limited to what the data can support: frequency-titration consistency,
   not personalized gamma entrainment validation.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
PROV = ROOT / "provenance"
WORKSPACE = Path(r"G:\Alzheimers Dementia - Digital Twin")
ETTER_FIG2N = (
    WORKSPACE
    / "real_experimental_data"
    / "etter_2019"
    / "41467_2019_13260_MOESM4_ESM_fig_2n.csv"
)
IACC_VC = (
    WORKSPACE
    / "real_experimental_data"
    / "iaccarino_2016"
    / "iaccarino_vc_cohort_with_endogenous_gamma.csv"
)
PHASE2_SENS = (
    WORKSPACE
    / "Journal of Neural Engineering"
    / "revised_resubmission"
    / "github_zenodo_deposit"
    / "data"
    / "phase2_weighting_sensitivity.json"
)

SEED = 42
RNG = np.random.default_rng(SEED)


def _utcnow() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_etter_fig2n(path: Path) -> pd.DataFrame:
    """Parse MOESM4 fig. 2n stim/baseline ratios exactly as published."""
    raw = pd.read_csv(path, header=None)
    # Locate header row containing 'Mouse'
    header_idx = None
    for i, row in raw.iterrows():
        vals = [str(v).strip().lower() for v in row.tolist()]
        if "mouse" in vals:
            header_idx = i
            break
    if header_idx is None:
        raise RuntimeError(f"Could not find Mouse header in {path}")

    records = []
    for i in range(header_idx + 1, len(raw)):
        row = raw.iloc[i].tolist()
        if pd.isna(row[0]):
            continue
        try:
            mouse_num = int(float(row[0]))
        except (TypeError, ValueError):
            continue

        def _f(x):
            if pd.isna(x) or str(x).strip() == "":
                return np.nan
            return float(x)

        r10, r40, r80 = _f(row[1]), _f(row[2]), _f(row[3])
        complete = bool(np.isfinite(r10) and np.isfinite(r40) and np.isfinite(r80))
        measured = {
            10.0: r10,
            40.0: r40,
            80.0: r80,
        }
        finite = {f: v for f, v in measured.items() if np.isfinite(v)}
        if finite:
            preferred = float(max(finite, key=finite.get))
        else:
            preferred = np.nan

        records.append(
            {
                "mouse_id": f"Etter_{mouse_num}",
                "source_mouse_label": mouse_num,
                "stim_baseline_ratio_10hz": r10,
                "stim_baseline_ratio_40hz": r40,
                "stim_baseline_ratio_80hz": r80,
                "complete_tri_frequency_10_40_80": complete,
                "preferred_stim_frequency_hz": preferred,
                "preferred_stim_frequency_definition": (
                    "argmax of measured stim/baseline ratios among available "
                    "frequencies in {10,40,80} Hz; NOT a continuous endogenous gamma peak"
                ),
                "distance_preferred_stim_from_40hz": (
                    abs(preferred - 40.0) if np.isfinite(preferred) else np.nan
                ),
                "source_article": "Etter et al. Nat Commun 2019",
                "source_doi": "10.1038/s41467-019-13260-9",
                "source_file": "41467_2019_13260_MOESM4_ESM.xlsx",
                "source_sheet_or_panel": "fig_2n (stim/baseline ratio)",
                "transcription_method": "direct CSV export from publisher MOESM4 spreadsheet",
                "group": "PVJ20+_ChETA",
            }
        )
    df = pd.DataFrame(records).sort_values("source_mouse_label").reset_index(drop=True)
    if len(df) != 10:
        raise RuntimeError(f"Expected 10 Etter mice, got {len(df)}")
    if int(df["complete_tri_frequency_10_40_80"].sum()) != 7:
        raise RuntimeError("Expected exactly 7 complete tri-frequency mice")
    return df


def build_iaccarino_table(path: Path) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Build Iaccarino VC ELISA table with raw Aβ and dark-baseline fractional reduction.

    Endogenous/intrinsic gamma columns are intentionally omitted from the animal table.
    """
    src = pd.read_csv(path)
    required = {
        "mouse_id",
        "stimulation_frequency",
        "amyloid_ab40_pg_ml",
        "treatment_response_40hz_raw",
        "experimental_group",
    }
    missing = required - set(src.columns)
    if missing:
        raise RuntimeError(f"Iaccarino source missing columns: {missing}")

    # Recover dark baseline mean from animals with positive fractional reduction.
    implied = []
    for _, r in src.iterrows():
        raw = float(r["treatment_response_40hz_raw"])
        ab = float(r["amyloid_ab40_pg_ml"])
        if 0.0 < raw < 1.0:
            implied.append(ab / (1.0 - raw))
    dark_mean = float(np.mean(implied)) if implied else float("nan")
    # Recompute fractional reduction for auditability
    frac = np.maximum(0.0, (dark_mean - src["amyloid_ab40_pg_ml"].astype(float)) / dark_mean)

    df = pd.DataFrame(
        {
            "mouse_id": src["mouse_id"],
            "experimental_group": src["experimental_group"],
            "stimulation_frequency_hz": src["stimulation_frequency"].astype(float),
            "stim_mismatch_from_40hz": (src["stimulation_frequency"].astype(float) - 40.0).abs(),
            "amyloid_ab40_pg_ml": src["amyloid_ab40_pg_ml"].astype(float),
            "dark_1hr_vc_baseline_mean_pg_ml": dark_mean,
            "fractional_ab40_reduction_vs_dark": frac,
            "no_reduction_vs_dark": frac <= 0.0,
            "group": "5xFAD",
            "source_article": "Iaccarino et al. Nature 2016",
            "source_doi": "10.1038/nature20587",
            "source_file": "NIHMS886326 (PMC Extended Data Table 1 / Aβ1-40 ELISA)",
            "source_sheet_or_panel": "Extended Data Table 1 — 20/40/80 Hz 1 hr VC arms",
            "transcription_method": "programmatic parse of NIHMS PDF Extended Data Table 1",
            "endogenous_gamma_note": (
                "NOT included: paired endogenous gamma was not published for these ELISA mice. "
                "Any SWR peaks from Ext. Data Fig. 1 are a separate electrophysiology cohort "
                "and must not be treated as per-animal predictors for this table."
            ),
        }
    )
    # Sanity: recomputed frac matches deposited raw reduction to numerical tolerance
    max_abs = float(np.max(np.abs(frac.values - src["treatment_response_40hz_raw"].astype(float).values)))
    meta = {
        "n": int(len(df)),
        "dark_1hr_vc_baseline_mean_pg_ml": dark_mean,
        "n_no_reduction_vs_dark": int((frac <= 0).sum()),
        "n_unique_ab40_values": int(df["amyloid_ab40_pg_ml"].nunique()),
        "recompute_max_abs_diff_vs_source_raw": max_abs,
        "arms": {
            str(int(k)): int(v)
            for k, v in df["stimulation_frequency_hz"].value_counts().sort_index().items()
        },
    }
    if max_abs > 1e-6:
        raise RuntimeError(f"Dark-baseline recompute mismatch: {max_abs}")
    if len(df) != 18:
        raise RuntimeError(f"Expected 18 Iaccarino animals, got {len(df)}")
    return df, meta


def spearman_with_permutation(
    x: np.ndarray, y: np.ndarray, n_perm: int = 100_000
) -> Dict[str, Any]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rho, p_asymp = stats.spearmanr(x, y)
    # Exact-ish two-sided permutation p-value
    count = 0
    abs_rho = abs(float(rho))
    for _ in range(n_perm):
        yp = RNG.permutation(y)
        r_perm, _ = stats.spearmanr(x, yp)
        if abs(float(r_perm)) >= abs_rho - 1e-15:
            count += 1
    p_perm = (count + 1) / (n_perm + 1)

    # Bootstrap CI for rho
    boots = []
    n = len(x)
    for _ in range(10_000):
        idx = RNG.integers(0, n, n)
        if len(np.unique(x[idx])) < 2 or len(np.unique(y[idx])) < 2:
            continue
        rb, _ = stats.spearmanr(x[idx], y[idx])
        if np.isfinite(rb):
            boots.append(float(rb))
    if boots:
        lo, hi = np.percentile(boots, [2.5, 97.5])
    else:
        lo = hi = float("nan")

    pr, pp = stats.pearsonr(x, y)
    return {
        "n": int(n),
        "spearman_rho": float(rho),
        "spearman_p_asymptotic_t": float(p_asymp),
        "spearman_p_permutation": float(p_perm),
        "spearman_rho_bootstrap_ci95": [float(lo), float(hi)],
        "pearson_r": float(pr),
        "pearson_p": float(pp),
        "n_permutations": n_perm,
        "n_bootstrap": len(boots),
        "random_seed": SEED,
    }


def power_spearman(rho_alt: float, n: int, alpha: float = 0.05) -> Dict[str, Any]:
    """Two-sided Fisher-z approximation (transparent; underpowered warning retained)."""
    if n < 4:
        return {"n": n, "approx_power": float("nan"), "note": "n too small"}
    z_r = np.arctanh(np.clip(rho_alt, -0.99, 0.99))
    z_crit = stats.norm.ppf(1 - alpha / 2)
    se = 1.0 / math.sqrt(n - 3)
    # Power = P(|Z| > z_crit) under N(z_r/se * se wait): noncentrality = z_r/se
    ncp = z_r / se
    power = float(stats.norm.cdf(-z_crit - ncp) + (1 - stats.norm.cdf(z_crit - ncp)))
    return {
        "n": n,
        "alpha": alpha,
        "rho_alternative": float(rho_alt),
        "approx_power_two_sided_fisher_z": power,
        "method": "Fisher z approximation applied to Spearman (approximate for small n)",
    }


def subgroup_contrast(df: pd.DataFrame, ycol: str) -> Dict[str, Any]:
    close = df.loc[df["stim_mismatch_from_40hz"] < 5, ycol].astype(float)
    far = df.loc[df["stim_mismatch_from_40hz"] >= 5, ycol].astype(float)
    # Mann-Whitney U (exact for small n where scipy supports)
    u, p = stats.mannwhitneyu(close, far, alternative="two-sided", method="auto")
    return {
        "close_n": int(len(close)),
        "far_n": int(len(far)),
        "close_mean": float(close.mean()),
        "far_mean": float(far.mean()),
        "close_median": float(close.median()),
        "far_median": float(far.median()),
        "mannwhitney_U": float(u),
        "mannwhitney_p": float(p),
        "definition": "close: |f_stim-40|<5 (40 Hz arm); far: |f_stim-40|>=5 (20 and 80 Hz arms)",
    }


def analysis_specification() -> Dict[str, Any]:
    """Locked analysis specification for the corrected reanalysis (NOT a backdated prereg)."""
    return {
        "document_type": "corrected_reanalysis_specification",
        "version": "2.0.0",
        "created_utc": _utcnow(),
        "honesty_notice": (
            "This is an analysis specification for the post-rejection corrected archive. "
            "It is NOT a prospective pre-registration fixed before seeing the data. "
            "The prior PHASE3_PREREGISTRATION.json incorrectly embedded finished coefficients "
            "and pass/fail verdicts inside the plan file; that practice is discontinued."
        ),
        "claim_boundary": {
            "supported": [
                "Frequency-titration consistency: in Iaccarino VC ELISA (n=18), greater assigned flicker mismatch from 40 Hz associates with lower Aβ1-40 reduction vs dark baseline.",
                "Etter complete tri-frequency subset (n=7) is underpowered; report point estimate with permutation CI; do not claim personalization validation.",
                "Population-level Brian2/digital-twin simulation is hypothesis-generating only.",
            ],
            "not_supported": [
                "Personalized gamma entrainment validated in individuals.",
                "Independent out-of-sample prediction of unseen experiments.",
                "Pooled n=28 ML as confirmatory evidence.",
                "28 Hz as a robust empirical stimulation optimum.",
                "Iaccarino per-animal endogenous gamma predictors.",
                "Etter plaque burden paired by mouse ID (values exist anonymously; ID link unverified).",
            ],
        },
        "endpoints": {
            "P1_primary": {
                "study": "Iaccarino_2016",
                "n": 18,
                "inclusion": "All VC 1 hr flicker animals in Extended Data Table 1 arms 20/40/80 Hz",
                "predictor": "abs(stimulation_frequency_hz - 40)",
                "outcome_primary": "fractional_ab40_reduction_vs_dark",
                "outcome_sensitivity": "amyloid_ab40_pg_ml (raw scale; higher=worse)",
                "tests": [
                    "Spearman rho with permutation p-value and bootstrap CI",
                    "Pearson r (descriptive)",
                    "Mann-Whitney close vs far on fractional reduction",
                ],
                "interpretation_limit": (
                    "Confirms published frequency titration (40 Hz superior to 20/80), "
                    "not an individualized frequency-mismatch mechanism."
                ),
            },
            "S1_secondary": {
                "study": "Etter_2019",
                "n": 7,
                "inclusion_rule_code": (
                    "stim_baseline_ratio_10hz.notna() & "
                    "stim_baseline_ratio_40hz.notna() & "
                    "stim_baseline_ratio_80hz.notna()"
                ),
                "excluded_explicitly": ["Etter_181", "Etter_183", "Etter_184"],
                "exclusion_reason": "Only 10 Hz stim/baseline reported; 40 and 80 Hz missing in MOESM4 fig_2n",
                "predictor": "abs(preferred_stim_frequency_hz - 40) where preferred=argmax of measured ratios",
                "outcome": "stim_baseline_ratio_40hz (raw, not min-max normalized)",
                "tests": [
                    "Spearman rho with permutation p-value and bootstrap CI",
                ],
                "power": "Compute Fisher-z approximate power at n=7 for |rho|=0.5 and 0.7 and for observed |rho|",
                "interpretation_limit": "Underpowered consistency check; not confirmatory of personalization.",
            },
            "excluded": [
                "pooled_n28_primary",
                "meta_learning_confirmatory",
                "etter_plaque_ID_aligned_amyloid_proxy",
                "iaccarino_endogenous_gamma_distance_primary",
                "incomplete_Etter_mice_with_10Hz_substituted_as_40Hz_response",
            ],
        },
        "phase2_composite_policy": {
            "report_both_scales": ["raw_component_scores", "zscored_across_sweep"],
            "usable_components_in_simplified_sweep": ["gc_total", "pac_modulation_index"],
            "zero_components_disclosed": ["gamma_restoration", "coherence_restoration"],
            "headline_optimum_rule": (
                "Primary reported optimum uses equal weights on usable (non-identically-zero) "
                "z-scored components. GC-only optima are sensitivity analyses only."
            ),
            "known_from_prior_sensitivity_json": {
                "equal_usable_optimum_hz": 42.0,
                "gc_only_optimum_hz": 28.0,
                "noise_bootstrap_fraction_optimum_28": 0.0,
            },
        },
    }


def run_analyses(etter: pd.DataFrame, iacc: pd.DataFrame) -> Dict[str, Any]:
    # P1
    p1_frac = spearman_with_permutation(
        iacc["stim_mismatch_from_40hz"].values,
        iacc["fractional_ab40_reduction_vs_dark"].values,
    )
    # On raw Aβ, mismatch should correlate POSITIVELY (higher Aβ when farther from 40)
    p1_raw = spearman_with_permutation(
        iacc["stim_mismatch_from_40hz"].values,
        iacc["amyloid_ab40_pg_ml"].values,
    )
    p1_contrast = subgroup_contrast(iacc, "fractional_ab40_reduction_vs_dark")

    # S1
    s1_df = etter.loc[etter["complete_tri_frequency_10_40_80"]].copy()
    if sorted(s1_df["mouse_id"].tolist()) != sorted(
        [
            "Etter_177",
            "Etter_236",
            "Etter_238",
            "Etter_255",
            "Etter_333",
            "Etter_334",
            "Etter_354",
        ]
    ):
        raise RuntimeError(f"Unexpected S1 membership: {s1_df['mouse_id'].tolist()}")

    s1 = spearman_with_permutation(
        s1_df["distance_preferred_stim_from_40hz"].values,
        s1_df["stim_baseline_ratio_40hz"].values,
    )
    s1_power = {
        "observed_abs_rho": abs(s1["spearman_rho"]),
        "power_at_observed_abs_rho": power_spearman(abs(s1["spearman_rho"]), 7),
        "power_at_medium_effect_rho_0.5": power_spearman(0.5, 7),
        "power_at_large_effect_rho_0.7": power_spearman(0.7, 7),
        "comment": (
            "n=7 yields far below 80% power even for large effects; "
            "non-significance is inconclusive, not evidence of no effect."
        ),
    }

    # Descriptive: incomplete mice are excluded from S1 (document their 10 Hz-only status)
    incomplete = etter.loc[~etter["complete_tri_frequency_10_40_80"], "mouse_id"].tolist()

    # Subset sensitivity for transparency (not a primary endpoint)
    near = s1_df.loc[s1_df["distance_preferred_stim_from_40hz"] == 0]
    # All far complete animals in S1
    far_complete = s1_df.loc[s1_df["distance_preferred_stim_from_40hz"] > 0]
    # Also evaluate all 3-near + 4-of-7-far including incomplete far mice for transparency
    all_near = etter.loc[etter["preferred_stim_frequency_hz"] == 40]
    all_far = etter.loc[etter["preferred_stim_frequency_hz"] != 40]
    # For sensitivity on S1 statistic structure: predictor only takes 0 or 30 in this dataset
    subset_note = {
        "s1_n": int(len(s1_df)),
        "s1_n_near_preferred_40": int(len(near)),
        "s1_n_far": int(len(far_complete)),
        "predictor_unique_values_in_s1": sorted(
            s1_df["distance_preferred_stim_from_40hz"].dropna().unique().tolist()
        ),
        "warning": (
            "With only two distinct predictor values (0 and 30 Hz), S1 is effectively a "
            "two-group rank comparison and cannot demonstrate a continuous individualized "
            "mismatch mechanism."
        ),
        "incomplete_mice_excluded_from_s1": incomplete,
        "n_all_far_including_incomplete": int(len(all_far)),
        "n_all_near": int(len(all_near)),
    }

    phase2 = {}
    if PHASE2_SENS.exists():
        phase2 = json.loads(PHASE2_SENS.read_text(encoding="utf-8"))

    return {
        "generated_utc": _utcnow(),
        "random_seed": SEED,
        "P1_iaccarino_frequency_titration": {
            "endpoint": "P1",
            "outcome_fractional_reduction_vs_dark": p1_frac,
            "outcome_raw_ab40_pg_ml": p1_raw,
            "close_vs_far": p1_contrast,
            "claim": (
                "Assigned flicker mismatch from 40 Hz associates with lower fractional "
                "Aβ1-40 reduction vs dark baseline (frequency titration consistency)."
            ),
            "not_a_claim": "Personalized endogenous-gamma matching.",
        },
        "S1_etter_complete_tri_frequency": {
            "endpoint": "S1",
            "included_mouse_ids": s1_df["mouse_id"].tolist(),
            "statistics_raw_40hz_ratio": s1,
            "power": s1_power,
            "design_notes": subset_note,
            "claim": "Underpowered secondary consistency check only.",
            "significant_at_0.05_permutation": bool(s1["spearman_p_permutation"] < 0.05),
        },
        "phase2_sensitivity_excerpt": {
            "source": str(PHASE2_SENS) if PHASE2_SENS.exists() else None,
            "equal_usable_optimum_hz": phase2.get("equal_usable", {}).get("optimum_frequency_hz"),
            "gc_only_optimum_hz": phase2.get("gc_only", {}).get("optimum_frequency_hz"),
            "noise_bootstrap_fraction_optimum_28": phase2.get("noise_floor_bootstrap", {}).get(
                "fraction_optimum_equals_28"
            ),
            "raw_score_28_vs_40": phase2.get("grid_clarification", {}).get("raw_score_28_vs_40"),
            "scientific_conclusion": (
                "Under equal weights on usable z-scored components the optimum is ~42 Hz. "
                "28 Hz appears only under GC-only weighting and is unstable under 1% noise "
                "(bootstrap fraction at 28 Hz = 0). Do not headline 28 Hz as a robust optimum."
            ),
        },
    }


def write_markdown_docs(
    etter: pd.DataFrame,
    iacc: pd.DataFrame,
    iacc_meta: Dict[str, Any],
    results: Dict[str, Any],
) -> None:
    correction = f"""# Correction and supersession notice (v2.0.0)

**Date (UTC):** {_utcnow()}  
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
"""
    (ROOT / "CORRECTION_NOTICE.md").write_text(correction, encoding="utf-8")

    provenance = f"""# Data provenance ledger (v2.0.0)

## Etter et al., 2019 (Nature Communications)

| Field | Value |
|------|------|
| Article | Etter et al. Nat Commun 2019 |
| DOI | 10.1038/s41467-019-13260-9 |
| Source file | 41467_2019_13260_MOESM4_ESM (fig_2n sheet/export) |
| Local source used | `{ETTER_FIG2N}` |
| n rows | {len(etter)} |
| Complete 10+40+80 | {int(etter['complete_tri_frequency_10_40_80'].sum())} |
| Incomplete (10 Hz only) | {', '.join(etter.loc[~etter['complete_tri_frequency_10_40_80'],'mouse_id'])} |
| Preferred frequency | argmax among measured stim/baseline ratios (proxy only) |
| Plaque burden | **Not deposited for confirmatory use** (anonymous % list; mouse-ID pairing unverified) |

## Iaccarino et al., 2016 (Nature)

| Field | Value |
|------|------|
| Article | Iaccarino et al. Nature 2016 |
| DOI | 10.1038/nature20587 |
| Source | NIHMS886326 Extended Data Table 1 (Aβ1-40 ELISA, VC 1 hr arms) |
| Local source used | `{IACC_VC}` |
| n rows | {len(iacc)} |
| Dark 1 hr VC baseline mean (pg/ml) | {iacc_meta['dark_1hr_vc_baseline_mean_pg_ml']:.6f} |
| Animals with no reduction vs dark | {iacc_meta['n_no_reduction_vs_dark']} |
| Unique raw Aβ values | {iacc_meta['n_unique_ab40_values']} |
| Endogenous gamma | **Excluded from animal table** (not published paired to ELISA mice) |

## Reproducibility

```bash
python code/build_verified_archive.py
python code/run_corrected_analyses.py
```

Random seed: {SEED}
"""
    (PROV / "DATA_PROVENANCE.md").write_text(provenance, encoding="utf-8")

    readme = f"""# ad-gamma-digital-twin — corrected verified archive v2.0.0

**Status:** Corrected scientific archive after JNE-111725 rejection.  
**Supersedes:** `vsanker21/ad-gamma-digital-twin-jne` release v1.0.0 / Zenodo 10.5281/zenodo.22746385.

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
"""
    (ROOT / "README.md").write_text(readme, encoding="utf-8")

    # Human-readable results summary
    p1 = results["P1_iaccarino_frequency_titration"]["outcome_fractional_reduction_vs_dark"]
    s1 = results["S1_etter_complete_tri_frequency"]["statistics_raw_40hz_ratio"]
    summary = f"""# Corrected analysis results summary

Generated: {results['generated_utc']}

## P1 (Iaccarino, n=18) — fractional Aβ reduction vs dark

- Spearman ρ = {p1['spearman_rho']:.6f}
- Permutation p = {p1['spearman_p_permutation']:.6g}
- Asymptotic p = {p1['spearman_p_asymptotic_t']:.6g}
- Bootstrap 95% CI = [{p1['spearman_rho_bootstrap_ci95'][0]:.3f}, {p1['spearman_rho_bootstrap_ci95'][1]:.3f}]
- Close vs far Mann–Whitney p = {results['P1_iaccarino_frequency_titration']['close_vs_far']['mannwhitney_p']:.6g}

## S1 (Etter complete tri-frequency, n=7) — raw 40 Hz ratio

- Included: {', '.join(results['S1_etter_complete_tri_frequency']['included_mouse_ids'])}
- Spearman ρ = {s1['spearman_rho']:.6f}
- Permutation p = {s1['spearman_p_permutation']:.6g}
- Bootstrap 95% CI = [{s1['spearman_rho_bootstrap_ci95'][0]:.3f}, {s1['spearman_rho_bootstrap_ci95'][1]:.3f}]
- Approx power at n=7, |ρ|=0.5: {results['S1_etter_complete_tri_frequency']['power']['power_at_medium_effect_rho_0.5']['approx_power_two_sided_fisher_z']:.3f}
- Approx power at n=7, |ρ|=0.7: {results['S1_etter_complete_tri_frequency']['power']['power_at_large_effect_rho_0.7']['approx_power_two_sided_fisher_z']:.3f}

## Phase 2

{results['phase2_sensitivity_excerpt']['scientific_conclusion']}
"""
    (RESULTS / "RESULTS_SUMMARY.md").write_text(summary, encoding="utf-8")


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    PROV.mkdir(parents=True, exist_ok=True)

    etter = parse_etter_fig2n(ETTER_FIG2N)
    iacc, iacc_meta = build_iaccarino_table(IACC_VC)

    etter_path = DATA / "etter_2019_moesm4_fig2n_raw.csv"
    iacc_path = DATA / "iaccarino_2016_ext_data_table1_raw.csv"
    etter.to_csv(etter_path, index=False)
    iacc.to_csv(iacc_path, index=False)

    spec = analysis_specification()
    (DATA / "ANALYSIS_SPECIFICATION.json").write_text(
        json.dumps(spec, indent=2), encoding="utf-8"
    )

    # Copy phase2 sensitivity if present
    if PHASE2_SENS.exists():
        (DATA / "phase2_weighting_sensitivity.json").write_text(
            PHASE2_SENS.read_text(encoding="utf-8"), encoding="utf-8"
        )

    results = run_analyses(etter, iacc)
    results["iaccarino_table_meta"] = iacc_meta
    (RESULTS / "CORRECTED_ANALYSIS_RESULTS.json").write_text(
        json.dumps(results, indent=2), encoding="utf-8"
    )

    # Verification fingerprint
    fingerprint = {
        "generated_utc": _utcnow(),
        "etter_sha256_n": len(etter),
        "etter_complete_n": int(etter["complete_tri_frequency_10_40_80"].sum()),
        "iaccarino_n": len(iacc),
        "iaccarino_unique_ab40": int(iacc["amyloid_ab40_pg_ml"].nunique()),
        "p1_rho": results["P1_iaccarino_frequency_titration"]["outcome_fractional_reduction_vs_dark"][
            "spearman_rho"
        ],
        "s1_rho": results["S1_etter_complete_tri_frequency"]["statistics_raw_40hz_ratio"][
            "spearman_rho"
        ],
        "s1_included": results["S1_etter_complete_tri_frequency"]["included_mouse_ids"],
    }
    (PROV / "verification_fingerprint.json").write_text(
        json.dumps(fingerprint, indent=2), encoding="utf-8"
    )

    write_markdown_docs(etter, iacc, iacc_meta, results)

    # Thin runner alias
    runner = '''"""Run corrected analyses from deposited raw tables (no rebuild from workspace)."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
SEED = 42
RNG = np.random.default_rng(SEED)

def spearman_perm(x, y, n_perm=100000):
    rho, p_asymp = stats.spearmanr(x, y)
    count = 0
    for _ in range(n_perm):
        r, _ = stats.spearmanr(x, RNG.permutation(y))
        if abs(r) >= abs(rho) - 1e-15:
            count += 1
    return float(rho), float(p_asymp), (count + 1) / (n_perm + 1)

def main():
    etter = pd.read_csv(ROOT / "data" / "etter_2019_moesm4_fig2n_raw.csv")
    iacc = pd.read_csv(ROOT / "data" / "iaccarino_2016_ext_data_table1_raw.csv")
    p1 = spearman_perm(iacc["stim_mismatch_from_40hz"], iacc["fractional_ab40_reduction_vs_dark"])
    s1 = etter[etter["complete_tri_frequency_10_40_80"]]
    s1s = spearman_perm(s1["distance_preferred_stim_from_40hz"], s1["stim_baseline_ratio_40hz"])
    out = {
        "P1_spearman_rho": p1[0],
        "P1_p_asymp": p1[1],
        "P1_p_perm": p1[2],
        "S1_n": int(len(s1)),
        "S1_mouse_ids": s1["mouse_id"].tolist(),
        "S1_spearman_rho": s1s[0],
        "S1_p_asymp": s1s[1],
        "S1_p_perm": s1s[2],
    }
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
'''
    (ROOT / "code" / "run_corrected_analyses.py").write_text(runner, encoding="utf-8")
    (ROOT / "code" / "requirements.txt").write_text(
        "numpy>=1.24\npandas>=2.0\nscipy>=1.10\n", encoding="utf-8"
    )
    (ROOT / "LICENSE").write_text(
        "MIT License\n\nCopyright (c) 2026 Vivek Sanker, Shriya Awasthi, Atman Desai, Sharon Sha\n\n"
        "Permission is hereby granted, free of charge, to any person obtaining a copy of this software "
        "and associated documentation files (the \"Software\"), to deal in the Software without restriction, "
        "including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, "
        "and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, "
        "subject to the following conditions:\n\n"
        "The above copyright notice and this permission notice shall be included in all copies or substantial "
        "portions of the Software.\n\n"
        "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED.\n",
        encoding="utf-8",
    )

    print(json.dumps(fingerprint, indent=2))


if __name__ == "__main__":
    main()
