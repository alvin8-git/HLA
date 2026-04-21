import pandas as pd
import numpy as np
import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))

from hwe_test import (
    compute_allele_frequencies,
    compute_hwe_stats,
    run_em_haplotypes,
    _enum_phase_configs,
    _em_full_phase,
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


# ---------------------------------------------------------------------------
# Tests for _enum_phase_configs
# ---------------------------------------------------------------------------

def test_enum_phase_configs_all_homozygous():
    geno = [('A', 'A'), ('B', 'B'), ('C', 'C'), ('D', 'D'), ('E', 'E')]
    configs = _enum_phase_configs(geno)
    assert len(configs) == 1
    h1, h2 = configs[0]
    assert h1 == h2 == ('A', 'B', 'C', 'D', 'E')


def test_enum_phase_configs_one_het_locus():
    # 1 het locus → 2^0 = 1 distinct phase
    geno = [('A', 'a'), ('B', 'B'), ('C', 'C'), ('D', 'D'), ('E', 'E')]
    configs = _enum_phase_configs(geno)
    assert len(configs) == 1
    h1, h2 = configs[0]
    assert set(h1) | set(h2)  # both haplotypes present
    assert h1[0] != h2[0]     # het locus differs


def test_enum_phase_configs_two_het_loci():
    # 2 het loci → 2^1 = 2 distinct phases (coupling and repulsion)
    geno = [('A', 'a'), ('B', 'b'), ('C', 'C'), ('D', 'D'), ('E', 'E')]
    configs = _enum_phase_configs(geno)
    assert len(configs) == 2
    hap_sets = [frozenset([h1, h2]) for h1, h2 in configs]
    assert hap_sets[0] != hap_sets[1], "Two phases should produce different diplotypes"


def test_enum_phase_configs_five_het_loci():
    # 5 het loci → 2^4 = 16 distinct phases
    geno = [('A', 'a'), ('B', 'b'), ('C', 'c'), ('D', 'd'), ('E', 'e')]
    configs = _enum_phase_configs(geno)
    assert len(configs) == 16
    # All configs are distinct
    assert len(set(frozenset([h1, h2]) for h1, h2 in configs)) == 16


def test_enum_phase_configs_covers_all_alleles():
    geno = [('X', 'Y'), ('P', 'Q'), ('A', 'A'), ('B', 'B'), ('C', 'C')]
    configs = _enum_phase_configs(geno)
    for h1, h2 in configs:
        for locus_idx, (a1, a2) in enumerate(geno):
            assert frozenset([h1[locus_idx], h2[locus_idx]]) == frozenset([a1, a2])


# ---------------------------------------------------------------------------
# Test that full EM recovers correct phase when column assignment is scrambled
# ---------------------------------------------------------------------------

def _make_phase_test_df(n_anchor=80, n_scrambled=40):
    """
    Construct a dataset with known true haplotypes h_true1 and h_true2.

    Anchoring individuals (n_anchor each): homozygous at all 5 loci for one
    true haplotype — no phase ambiguity, establishes LD signal in the EM.

    Scrambled individuals (n_scrambled): heterozygous at all 5 loci (true
    phase is h_true1/h_true2) but the HLA-B allele column is swapped,
    creating a mis-stated phase in the raw data.

    The proper EM uses the LD established by the anchors to correctly resolve
    the scrambled individuals back to h_true1 and h_true2.
    """
    LOCI_NAMES = ['HLA-A', 'HLA-B', 'HLA-C', 'DRB1', 'DQB1']
    h_true1 = ('02:01', '46:01', '01:02', '09:01', '03:03')
    h_true2 = ('11:01', '40:01', '03:02', '03:01', '02:01')

    rows = []
    sid = 0

    def add_individual(allele_pairs):
        nonlocal sid
        for locus, (a1, a2) in zip(LOCI_NAMES, allele_pairs):
            rows.append({'sample_id': f'S{sid}', 'source': 'BMDP_OUT',
                         'ethnicity': 'Chinese', 'locus': locus,
                         'allele1': a1, 'allele2': a2})
        sid += 1

    # Anchors: homozygous for h_true1 and h_true2 (unambiguous phase)
    for _ in range(n_anchor):
        add_individual([(a, a) for a in h_true1])
    for _ in range(n_anchor):
        add_individual([(a, a) for a in h_true2])

    # Scrambled: all loci het (h_true1/h_true2), but HLA-B column is swapped
    for _ in range(n_scrambled):
        pairs = list(zip(h_true1, h_true2))
        pairs[1] = (h_true2[1], h_true1[1])   # swap HLA-B col assignment
        add_individual(pairs)

    return pd.DataFrame(rows), h_true1, h_true2


def test_full_em_resolves_phase_correctly():
    """
    The proper EM must identify the two true haplotypes as dominant even when
    some individuals have a scrambled (wrong-phase) column assignment.

    The anchoring homozygous individuals establish LD so the E-step assigns
    near-zero weight to the spurious 'scrambled' phase configurations.
    """
    from hwe_test import _pivot_to_individuals
    df, h_true1, h_true2 = _make_phase_test_df()
    wide = _pivot_to_individuals(df[df['ethnicity'] == 'Chinese'])
    result = _em_full_phase(wide, max_iter=300, tol=1e-8)

    sorted_haps = sorted(result.items(), key=lambda x: -x[1])
    top2 = {h for h, _ in sorted_haps[:2]}
    assert h_true1 in top2, f"h_true1 not in top-2: {sorted_haps[:5]}"
    assert h_true2 in top2, f"h_true2 not in top-2: {sorted_haps[:5]}"
    assert result[h_true1] > 0.45
    assert result[h_true2] > 0.45
