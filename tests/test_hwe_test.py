import pandas as pd
import numpy as np
import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))

from hwe_test import (
    compute_allele_frequencies,
    compute_hwe_stats,
    run_em_haplotypes,
)

def _make_typed_df(n=100, seed=0):
    """Make a small tidy HLA df with 5 loci, all typed."""
    rng = np.random.default_rng(seed)
    alleles_a = ['01:01', '02:01', '03:01', '24:02']
    alleles_b = ['07:02', '08:01', '35:01', '57:01']
    alleles_c = ['03:04', '04:01', '06:02', '07:01']
    alleles_d = ['01:01', '03:01', '04:01', '15:01']
    alleles_q = ['02:01', '03:01', '05:01', '06:01']
    loci_alleles = {
        'HLA-A': alleles_a, 'HLA-B': alleles_b, 'HLA-C': alleles_c,
        'DRB1': alleles_d, 'DQB1': alleles_q,
    }
    rows = []
    for sid in range(n):
        for locus, opts in loci_alleles.items():
            a1 = rng.choice(opts)
            a2 = rng.choice(opts)
            rows.append({'sample_id': f'S{sid}', 'source': 'BMDP_OUT',
                         'ethnicity': 'Chinese', 'locus': locus,
                         'allele1': a1, 'allele2': a2})
    return pd.DataFrame(rows)


def test_compute_allele_frequencies_sums_to_one():
    df = _make_typed_df(100)
    result = compute_allele_frequencies(df)
    # Each (ethnicity, locus) group should sum to ~1.0
    grouped = result.groupby(['ethnicity', 'locus'])['frequency'].sum()
    for val in grouped:
        assert abs(val - 1.0) < 1e-9, f"Frequencies don't sum to 1: {val}"


def test_compute_allele_frequencies_uses_both_alleles():
    """One sample, HLA-A: allele1=01:01, allele2=02:01 → each gets 0.5."""
    df = pd.DataFrame([{
        'sample_id': 'S1', 'source': 'BMDP_OUT', 'ethnicity': 'Chinese',
        'locus': 'HLA-A', 'allele1': '01:01', 'allele2': '02:01'
    }])
    result = compute_allele_frequencies(df)
    freqs = result.set_index('allele')['frequency'].to_dict()
    assert abs(freqs['01:01'] - 0.5) < 1e-9
    assert abs(freqs['02:01'] - 0.5) < 1e-9


def test_compute_allele_frequencies_ignores_nan_allele2():
    """allele2=NaN should not affect frequency calculation."""
    df = pd.DataFrame([
        {'sample_id': 'S1', 'source': 'BMDP_OUT', 'ethnicity': 'Chinese',
         'locus': 'HLA-A', 'allele1': '01:01', 'allele2': np.nan},
        {'sample_id': 'S2', 'source': 'BMDP_OUT', 'ethnicity': 'Chinese',
         'locus': 'HLA-A', 'allele1': '01:01', 'allele2': '02:01'},
    ])
    result = compute_allele_frequencies(df)
    freqs = result.set_index('allele')['frequency'].to_dict()
    # 3 allele observations: 01:01 appears twice (S1 allele1 + S2 allele1), 02:01 once
    assert abs(freqs['01:01'] - 2/3) < 1e-9
    assert abs(freqs['02:01'] - 1/3) < 1e-9


def test_compute_hwe_stats_structure():
    df = _make_typed_df(200)
    allele_freqs = compute_allele_frequencies(df)
    result = compute_hwe_stats(df, allele_freqs)
    assert set(result.columns) >= {'ethnicity', 'locus', 'n_individuals',
                                    'H_obs', 'H_exp', 'chi2_stat', 'p_value', 'significant'}
    assert len(result) == 5  # 5 loci for 1 ethnicity
    assert (result['H_exp'] >= 0).all() and (result['H_exp'] <= 1).all()
    assert (result['H_obs'] >= 0).all() and (result['H_obs'] <= 1).all()


def test_run_em_haplotypes_returns_valid_frequencies():
    df = _make_typed_df(50)
    result = run_em_haplotypes(df)
    assert 'ethnicity' in result.columns
    assert 'haplotype' in result.columns
    assert 'frequency' in result.columns
    # All frequencies positive and sum <= 1 per ethnicity (may not sum to 1 after threshold)
    assert (result['frequency'] > 0).all()
    total = result.groupby('ethnicity')['frequency'].sum().iloc[0]
    assert total <= 1.0 + 1e-6
    # Haplotype format: 5 alleles pipe-separated
    sample_hap = result['haplotype'].iloc[0]
    assert sample_hap.count('|') == 4
