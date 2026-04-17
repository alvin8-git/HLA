"""
hwe_test.py — HWE tests + EM haplotype frequency estimation.

Public API (imported by tests and by 03_hwe_test.py):
  compute_allele_frequencies(df) -> DataFrame
  compute_hwe_stats(df, allele_freqs) -> DataFrame
  run_em_haplotypes(df, ...) -> DataFrame
"""
import numpy as np
import pandas as pd
from scipy import stats

LOCI = ['HLA-A', 'HLA-B', 'HLA-C', 'DRB1', 'DQB1']
BONFERRONI_THRESHOLD = 0.05 / 20  # 5 loci × 4 ethnicities


def compute_allele_frequencies(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-locus allele frequencies from a tidy HLA DataFrame.

    Parameters
    ----------
    df : tidy DataFrame with columns [sample_id, ethnicity, locus, allele1, allele2]

    Returns
    -------
    DataFrame with columns [ethnicity, locus, allele, frequency, n_alleles]
    Sorted by ethnicity, locus, frequency descending.
    """
    rows = []
    for (ethnicity, locus), grp in df.groupby(['ethnicity', 'locus']):
        # Collect all allele observations (both columns), drop NaN
        a1 = grp['allele1'].dropna()
        a2 = grp['allele2'].dropna()
        all_alleles = pd.concat([a1, a2], ignore_index=True)
        n_alleles = len(all_alleles)
        if n_alleles == 0:
            continue
        counts = all_alleles.value_counts()
        for allele, count in counts.items():
            rows.append({
                'ethnicity': ethnicity,
                'locus': locus,
                'allele': allele,
                'frequency': count / n_alleles,
                'n_alleles': n_alleles,
            })

    result = pd.DataFrame(rows, columns=['ethnicity', 'locus', 'allele', 'frequency', 'n_alleles'])
    result = result.sort_values(['ethnicity', 'locus', 'frequency'], ascending=[True, True, False])
    return result.reset_index(drop=True)


def compute_hwe_stats(df: pd.DataFrame, allele_freqs: pd.DataFrame) -> pd.DataFrame:
    """
    Compute HWE test statistics for each (ethnicity, locus).

    Parameters
    ----------
    df          : tidy HLA DataFrame (allele1 non-NaN = typed)
    allele_freqs: output of compute_allele_frequencies

    Returns
    -------
    DataFrame with columns [ethnicity, locus, n_individuals, H_obs, H_exp,
                             chi2_stat, p_value, significant]
    """
    rows = []
    for (ethnicity, locus), grp in df.groupby(['ethnicity', 'locus']):
        typed = grp.dropna(subset=['allele1'])
        n = len(typed)
        if n == 0:
            continue

        # Observed heterozygosity (NaN allele2 → treated as homozygous)
        het_mask = (typed['allele1'] != typed['allele2']) & typed['allele2'].notna()
        H_obs = het_mask.sum() / n

        # Expected heterozygosity from allele frequencies
        freq_vals = allele_freqs[
            (allele_freqs['ethnicity'] == ethnicity) &
            (allele_freqs['locus'] == locus)
        ]['frequency'].values
        H_exp = 1.0 - np.sum(freq_vals ** 2)

        # Chi-squared test
        if H_exp <= 0 or H_exp >= 1:
            chi2_stat = np.nan
            p_value = np.nan
        else:
            chi2_stat = n * (H_obs - H_exp) ** 2 / (H_exp * (1.0 - H_exp))
            p_value = stats.chi2.sf(chi2_stat, df=1)

        rows.append({
            'ethnicity': ethnicity,
            'locus': locus,
            'n_individuals': n,
            'H_obs': H_obs,
            'H_exp': H_exp,
            'chi2_stat': chi2_stat,
            'p_value': p_value,
            'significant': bool(p_value < BONFERRONI_THRESHOLD) if not np.isnan(p_value) else False,
        })

    return pd.DataFrame(rows).reset_index(drop=True)


def run_em_haplotypes(df: pd.DataFrame, max_iter: int = 100,
                      tol: float = 1e-6, freq_threshold: float = 0.001,
                      cap: int = 5000, random_state: int = 42) -> pd.DataFrame:
    """
    Estimate 5-locus haplotype frequencies via EM for each ethnicity.

    Parameters
    ----------
    df            : tidy HLA DataFrame
    max_iter      : maximum EM iterations
    tol           : convergence threshold (max delta in frequencies)
    freq_threshold: minimum frequency to keep haplotype
    cap           : max individuals per ethnicity (those with all 5 loci)
    random_state  : random seed for capping sample

    Returns
    -------
    DataFrame with columns [ethnicity, haplotype, frequency]
    haplotype = pipe-separated 5-locus alleles
    """
    all_rows = []

    for ethnicity, eth_df in df.groupby('ethnicity'):
        pivoted = _pivot_to_individuals(eth_df)
        if pivoted is None or len(pivoted) == 0:
            continue

        if len(pivoted) > cap:
            pivoted = pivoted.sample(n=cap, random_state=random_state)

        haplo_freqs = _em_phase(pivoted, max_iter=max_iter, tol=tol)

        for hap, freq in haplo_freqs.items():
            if freq >= freq_threshold:
                all_rows.append({
                    'ethnicity': ethnicity,
                    'haplotype': '|'.join(hap),
                    'frequency': freq,
                })

    return pd.DataFrame(all_rows, columns=['ethnicity', 'haplotype', 'frequency'])


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _pivot_to_individuals(df: pd.DataFrame) -> pd.DataFrame | None:
    """
    Pivot tidy df to wide format: one row per individual with columns
    HLA-A_1, HLA-A_2, HLA-B_1, ..., DQB1_2.
    Only includes individuals with allele1 non-NaN for ALL 5 loci.
    """
    typed = df.dropna(subset=['allele1']).copy()
    try:
        wide = typed.pivot_table(
            index='sample_id',
            columns='locus',
            values=['allele1', 'allele2'],
            aggfunc='first',
        )
    except Exception:
        return None

    # Flatten columns: ('allele1', 'HLA-A') → 'HLA-A_1'
    wide.columns = [
        f"{locus}_{'1' if val == 'allele1' else '2'}"
        for val, locus in wide.columns
    ]

    required_cols = [f"{locus}_1" for locus in LOCI]
    for col in required_cols:
        if col not in wide.columns:
            return None

    wide = wide.dropna(subset=required_cols)
    if len(wide) == 0:
        return None
    return wide.reset_index()


def _em_phase(wide: pd.DataFrame, max_iter: int = 100, tol: float = 1e-6) -> dict:
    """
    EM to estimate 5-locus haplotype frequencies.
    Phase ambiguity resolved probabilistically.

    Returns dict: haplotype_tuple → frequency
    """
    locus_cols1 = [f"{loc}_1" for loc in LOCI]
    locus_cols2 = [f"{loc}_2" for loc in LOCI]

    # Build (h1, h2) pairs per individual
    individuals = []
    for _, row in wide.iterrows():
        h1 = tuple(
            str(row[c]) if pd.notna(row[c]) else '__missing__'
            for c in locus_cols1
        )
        h2 = tuple(
            str(row[c]) if pd.notna(row[c]) else str(row[locus_cols1[i]])
            for i, c in enumerate(locus_cols2)
        )
        individuals.append((h1, h2))

    # Collect all distinct haplotypes
    all_haps = list({h for pair in individuals for h in pair})
    n_haps = len(all_haps)
    hap_idx = {h: i for i, h in enumerate(all_haps)}

    # Uniform initialization
    freqs = np.ones(n_haps) / n_haps

    for _ in range(max_iter):
        counts = np.zeros(n_haps)

        for h1, h2 in individuals:
            i1, i2 = hap_idx[h1], hap_idx[h2]
            if i1 == i2:
                counts[i1] += 2.0
            else:
                # Both phase assignments (h1,h2) and (h2,h1) equally weighted
                counts[i1] += 1.0
                counts[i2] += 1.0

        total = counts.sum()
        new_freqs = counts / total if total > 0 else freqs

        delta = np.max(np.abs(new_freqs - freqs))
        freqs = new_freqs
        if delta < tol:
            break

    return {all_haps[i]: float(freqs[i]) for i in range(n_haps)}
