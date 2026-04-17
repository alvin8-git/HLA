"""
03_hwe_test.py — Step 3 of HLA pipeline.
Runs HWE tests and EM haplotype estimation, saving 3 output CSVs.
"""
import os
import sys
import pandas as pd

# Ensure analysis/ is on sys.path so hwe_test is importable
sys.path.insert(0, os.path.dirname(__file__))

from hwe_test import (
    compute_allele_frequencies,
    compute_hwe_stats,
    run_em_haplotypes,
    BONFERRONI_THRESHOLD,
)


def main():
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    input_path = os.path.join(data_dir, 'hla_clean.csv')

    print(f"Loading {input_path} ...")
    df_raw = pd.read_csv(input_path)

    # Filter to BMDP_OUT and SCBB_OUT sources only
    df = df_raw[df_raw['source'].isin(['BMDP_OUT', 'SCBB_OUT'])].copy()

    print(f"Filtered dataset: {len(df):,} rows, {df['sample_id'].nunique():,} unique samples")
    print(f"Ethnicities: {sorted(df['ethnicity'].unique())}")
    print(f"Loci: {sorted(df['locus'].unique())}")

    # Drop rows where allele1 is NaN for frequency/HWE computations
    df_typed = df.dropna(subset=['allele1']).copy()

    # 1. Allele frequencies
    print("\nComputing allele frequencies ...")
    allele_freqs = compute_allele_frequencies(df_typed)
    allele_freq_path = os.path.join(data_dir, 'allele_freqs_per_locus.csv')
    allele_freqs.to_csv(allele_freq_path, index=False)
    print(f"  Saved {len(allele_freqs):,} rows → {allele_freq_path}")

    # 2. HWE statistics
    print("\nComputing HWE statistics ...")
    hwe_results = compute_hwe_stats(df_typed, allele_freqs)
    hwe_path = os.path.join(data_dir, 'hwe_results.csv')
    hwe_results.to_csv(hwe_path, index=False)
    print(f"  Saved {len(hwe_results):,} rows → {hwe_path}")

    # 3. EM haplotype frequencies (pass full df including NaN allele2 rows —
    #    _pivot_to_individuals inside drops individuals missing any locus)
    print("\nRunning EM haplotype estimation (5-locus, per ethnicity) ...")
    haplo_freqs = run_em_haplotypes(df)
    haplo_path = os.path.join(data_dir, 'haplo_freqs_em.csv')
    haplo_freqs.to_csv(haplo_path, index=False)
    print(f"  Saved {len(haplo_freqs):,} haplotypes → {haplo_path}")

    # Summary
    print("\n=== SUMMARY ===")
    print(f"{'Ethnicity':<12} {'N_haplotypes':>14} {'HWE_violations':>16}")
    print("-" * 44)
    for ethnicity in sorted(df['ethnicity'].unique()):
        n_hap = len(haplo_freqs[haplo_freqs['ethnicity'] == ethnicity])
        hwe_eth = hwe_results[hwe_results['ethnicity'] == ethnicity]
        n_viol = int(hwe_eth['significant'].sum()) if len(hwe_eth) > 0 else 0
        print(f"{ethnicity:<12} {n_hap:>14,} {n_viol:>16}")

    print(f"\nBonferroni threshold: p < {BONFERRONI_THRESHOLD:.4f} (0.05 / 20)")
    sig = hwe_results[hwe_results['significant']]
    if len(sig) > 0:
        print(f"\nHWE violations ({len(sig)}):")
        for _, row in sig.iterrows():
            print(f"  {row['ethnicity']} / {row['locus']}: "
                  f"H_obs={row['H_obs']:.4f}, H_exp={row['H_exp']:.4f}, "
                  f"p={row['p_value']:.2e}")
    else:
        print("No HWE violations detected.")


if __name__ == '__main__':
    main()
