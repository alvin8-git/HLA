import pandas as pd
import numpy as np
import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))

from registry_model import (
    collapse_to_4locus,
    get_diplotype_frequencies,
    compute_coverage,
    find_registry_size,
    get_combined_haplotype_freqs,
    SG_WEIGHTS,
)


def _haplo_df(haplotypes_freqs):
    """Helper: build a minimal haplotype DataFrame."""
    rows = [{'ethnicity': 'Chinese', 'haplotype': h, 'frequency': f}
            for h, f in haplotypes_freqs]
    return pd.DataFrame(rows)


def test_collapse_to_4locus_strips_dqb1():
    df = _haplo_df([
        ('01:01|07:02|03:04|01:01|02:01', 0.6),
        ('01:01|07:02|03:04|01:01|03:01', 0.4),
    ])
    result = collapse_to_4locus(df)
    assert set(result['haplotype']) == {'01:01|07:02|03:04|01:01'}
    row = result[result['haplotype'] == '01:01|07:02|03:04|01:01'].iloc[0]
    assert abs(row['frequency'] - 1.0) < 1e-9


def test_get_diplotype_frequencies_sums_to_one():
    df = _haplo_df([('A|B|C|D|Q', 0.6), ('A2|B2|C2|D2|Q2', 0.4)])
    result = get_diplotype_frequencies(df)
    assert abs(result['diplotype_freq'].sum() - 1.0) < 1e-6


def test_get_diplotype_frequencies_homozygous():
    """Single haplotype with freq=1.0 → one diplotype with freq=1.0."""
    df = _haplo_df([('A|B|C|D|Q', 1.0)])
    result = get_diplotype_frequencies(df)
    assert len(result) == 1
    assert abs(result.iloc[0]['diplotype_freq'] - 1.0) < 1e-9
    assert result.iloc[0]['haplotype1'] == result.iloc[0]['haplotype2']


def test_compute_coverage_monotone_increasing():
    """Coverage should increase with N."""
    df = _haplo_df([('A|B|C|D|Q', 0.3), ('A2|B2|C2|D2|Q2', 0.7)])
    diplo = get_diplotype_frequencies(df)
    coverages = [compute_coverage(diplo, n) for n in [100, 1000, 10000]]
    assert coverages[0] < coverages[1] < coverages[2]
    assert all(0 <= c <= 1 for c in coverages)


def test_find_registry_size_reaches_target():
    """find_registry_size should return N where coverage >= target."""
    df = _haplo_df([('A|B|C|D|Q', 0.3), ('A2|B2|C2|D2|Q2', 0.7)])
    diplo = get_diplotype_frequencies(df)
    for target in [0.75, 0.85, 0.90]:
        n = find_registry_size(diplo, target)
        cov = compute_coverage(diplo, n)
        assert cov >= target - 1e-4, f"coverage {cov} < target {target} at N={n}"
        # verify N-1 is below target (we found minimum)
        if n > 1:
            cov_minus1 = compute_coverage(diplo, n - 1)
            assert cov_minus1 < target + 1e-4


def test_get_combined_haplotype_freqs_sums_to_one():
    """Combined frequencies should sum to 1.0."""
    dfs = {
        'Chinese': pd.DataFrame([{'haplotype': 'H1', 'frequency': 1.0}]),
        'Malay':   pd.DataFrame([{'haplotype': 'H1', 'frequency': 1.0}]),
        'Indian':  pd.DataFrame([{'haplotype': 'H1', 'frequency': 1.0}]),
        'Others':  pd.DataFrame([{'haplotype': 'H1', 'frequency': 1.0}]),
    }
    result = get_combined_haplotype_freqs(dfs)
    assert abs(result['frequency'].sum() - 1.0) < 1e-9
