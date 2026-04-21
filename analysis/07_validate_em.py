"""
07_validate_em.py
Validate our full-EM haplotype frequencies against Gene[Rate] estimates
from BMDPnSCBB.results.xlsx (Haplotype.* sheets).

Outputs:
  analysis/data/em_validation.csv   — matched haplotype comparison per ethnicity
  analysis/data/em_validation_summary.csv — per-ethnicity Spearman r, RMSE
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA_DIR = os.path.join(HERE, 'data')
XLSX = os.path.join(ROOT, 'BMDPnSCBB.results.xlsx')

ETHNICITIES = ['Chinese', 'Malay', 'Indian', 'Others']

# Gene[Rate] haplotype format: A*33:03~B*58:01~C*03:02~DRB1*03:01~DQB1*02:01
# Our EM format:                  33:03|58:01|03:02|03:01|02:01
LOCUS_PREFIXES = ['A*', 'B*', 'C*', 'DRB1*', 'DQB1*']


def parse_genrate_sheet(eth: str) -> pd.DataFrame:
    """Parse a Gene[Rate] Haplotype.* sheet → DataFrame with [haplotype, freq_genrate]."""
    df = pd.read_excel(XLSX, sheet_name=f'Haplotype.{eth}', header=None)
    # Row 0 is header; col 0 = haplotype string, col 1 = "Est %"
    df = df.iloc[1:].copy()
    df.columns = ['haplotype_raw'] + [f'col{i}' for i in range(1, len(df.columns))]
    df = df[['haplotype_raw', 'col1']].rename(columns={'col1': 'freq_genrate'})
    df = df.dropna(subset=['haplotype_raw', 'freq_genrate'])
    df['freq_genrate'] = pd.to_numeric(df['freq_genrate'], errors='coerce')
    df = df.dropna(subset=['freq_genrate'])

    # Convert Gene[Rate] format to our pipe-separated format
    def convert(s):
        parts = str(s).split('~')
        if len(parts) != 5:
            return None
        alleles = []
        for part, prefix in zip(parts, LOCUS_PREFIXES):
            if part.startswith(prefix):
                alleles.append(part[len(prefix):])
            else:
                alleles.append(part)
        return '|'.join(alleles)

    df['haplotype'] = df['haplotype_raw'].apply(convert)
    df = df.dropna(subset=['haplotype'])
    return df[['haplotype', 'freq_genrate']].reset_index(drop=True)


def main():
    em_df = pd.read_csv(os.path.join(DATA_DIR, 'haplo_freqs_em.csv'))

    all_matched = []
    summary_rows = []

    for eth in ETHNICITIES:
        print(f"\n=== {eth} ===")
        gr = parse_genrate_sheet(eth)
        em = em_df[em_df['ethnicity'] == eth][['haplotype', 'frequency']].rename(
            columns={'frequency': 'freq_em'})

        merged = pd.merge(gr, em, on='haplotype', how='outer')
        merged['ethnicity'] = eth
        n_gr = gr['haplotype'].nunique()
        n_em = em['haplotype'].nunique()
        n_both = merged.dropna(subset=['freq_genrate', 'freq_em']).shape[0]
        n_gr_only = merged[merged['freq_em'].isna()].shape[0]
        n_em_only = merged[merged['freq_genrate'].isna()].shape[0]

        print(f"  Gene[Rate] haplotypes : {n_gr}")
        print(f"  Our EM haplotypes     : {n_em}")
        print(f"  Matched (both)        : {n_both}")
        print(f"  Gene[Rate] only       : {n_gr_only}")
        print(f"  Our EM only           : {n_em_only}")

        matched = merged.dropna(subset=['freq_genrate', 'freq_em']).copy()
        if len(matched) >= 5:
            r, p = spearmanr(matched['freq_genrate'], matched['freq_em'])
            rmse = np.sqrt(np.mean((matched['freq_genrate'] - matched['freq_em']) ** 2))
            freq_matched_gr = matched['freq_genrate'].sum()
            freq_matched_em = matched['freq_em'].sum()
            print(f"  Spearman r            : {r:.4f}  (p={p:.2e})")
            print(f"  RMSE                  : {rmse:.5f}")
            print(f"  Coverage matched (GR) : {freq_matched_gr:.3f}")
            print(f"  Coverage matched (EM) : {freq_matched_em:.3f}")

            top10 = matched.sort_values('freq_genrate', ascending=False).head(10)
            print(f"\n  Top-10 Gene[Rate] haplotypes vs our EM:")
            print(f"  {'Haplotype':<45} {'GR':>7} {'EM':>7} {'delta':>8}")
            for _, row in top10.iterrows():
                delta = row['freq_em'] - row['freq_genrate']
                print(f"  {row['haplotype']:<45} {row['freq_genrate']:>7.4f} {row['freq_em']:>7.4f} {delta:>+8.4f}")

            summary_rows.append({
                'ethnicity': eth,
                'n_genrate': n_gr,
                'n_em': n_em,
                'n_matched': n_both,
                'n_gr_only': n_gr_only,
                'n_em_only': n_em_only,
                'spearman_r': round(r, 4),
                'spearman_p': f'{p:.2e}',
                'rmse': round(rmse, 6),
                'freq_coverage_gr': round(freq_matched_gr, 4),
                'freq_coverage_em': round(freq_matched_em, 4),
            })
        else:
            print("  Too few matched haplotypes for correlation.")

        all_matched.append(merged)

    all_df = pd.concat(all_matched, ignore_index=True)
    all_df.to_csv(os.path.join(DATA_DIR, 'em_validation.csv'), index=False)
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(DATA_DIR, 'em_validation_summary.csv'), index=False)
    print(f"\nSaved: {os.path.join(DATA_DIR, 'em_validation.csv')}")
    print(f"Saved: {os.path.join(DATA_DIR, 'em_validation_summary.csv')}")
    print("\n=== Summary ===")
    print(summary_df.to_string(index=False))


if __name__ == '__main__':
    main()
