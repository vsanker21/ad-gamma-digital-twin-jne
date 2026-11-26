"""Run corrected analyses from deposited raw tables (no rebuild from workspace)."""
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
