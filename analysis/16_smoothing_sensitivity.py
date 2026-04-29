"""
16_smoothing_sensitivity.py
Rare-haplotype smoothing sensitivity analysis.
Applies Laplace (add-α) pseudocount smoothing to EM haplotype frequencies
and compares N* with and without smoothing for all CMIO groups.

Outputs:
  data/smoothing_sensitivity.csv
"""
import os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from registry_model import get_diplotype_frequencies, find_registry_size

DATA_DIR = os.path.join(HERE, 'data')

THRESHOLDS  = [0.75, 0.85, 0.90, 0.95]
PSEUDOCOUNT = 0.001   # Dirichlet prior α added to each observed haplotype
ETHNICITIES = ['Chinese', 'Malay', 'Indian', 'Others']


def smooth_haplotypes(haplo_df: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """
    Add pseudocount alpha to each haplotype frequency, then renormalise.
    This assigns non-zero probability to all observed haplotypes,
    shrinking their frequencies toward uniformity and raising the effective
    probability mass available for rare diplotypes.
    """
    df = haplo_df.copy()
    df['frequency'] = df['frequency'] + alpha
    df['frequency'] /= df['frequency'].sum()
    return df


def main():
    haplo_df = pd.read_csv(os.path.join(DATA_DIR, 'haplo_freqs_em.csv'))
    rows = []

    for eth in ETHNICITIES:
        eth_df = haplo_df[haplo_df['ethnicity'] == eth].copy()
        if eth_df.empty:
            continue

        # Normalise original frequencies
        eth_df['frequency'] = eth_df['frequency'] / eth_df['frequency'].sum()
        eth_smooth = smooth_haplotypes(eth_df, PSEUDOCOUNT)

        diplo_orig   = get_diplotype_frequencies(eth_df)
        diplo_smooth = get_diplotype_frequencies(eth_smooth)

        for thr in THRESHOLDS:
            n_orig   = find_registry_size(diplo_orig,   thr)
            n_smooth = find_registry_size(diplo_smooth, thr)
            pct_diff = (n_smooth - n_orig) / n_orig * 100

            rows.append({
                'ethnicity':        eth,
                'target_coverage':  thr,
                'n_star_original':  n_orig,
                'n_star_smoothed':  n_smooth,
                'pct_increase':     round(pct_diff, 2),
            })
            print(f"  {eth:8s}  {thr:.0%}  orig={n_orig:,}  smooth={n_smooth:,}  "
                  f"Δ={pct_diff:+.2f}%")

    df = pd.DataFrame(rows)
    out = os.path.join(DATA_DIR, 'smoothing_sensitivity.csv')
    df.to_csv(out, index=False)
    print(f"\nSaved: {out}")

    print("\n=== Summary at 95% coverage ===")
    mask95 = df['target_coverage'] == 0.95
    print(df[mask95][['ethnicity', 'n_star_original', 'n_star_smoothed', 'pct_increase']]
          .to_string(index=False))


if __name__ == '__main__':
    main()
