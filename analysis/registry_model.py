"""
registry_model.py — Core math functions for HLA registry coverage modeling.
"""

import numpy as np
import pandas as pd
import math

# Singapore population ethnic weights
SG_WEIGHTS = {'Chinese': 0.77, 'Malay': 0.08, 'Indian': 0.09, 'Others': 0.06}


def collapse_to_4locus(haplotype_df):
    """
    Collapse 5-locus haplotypes (A|B|C|DRB1|DQB1) to 4-locus (A|B|C|DRB1)
    by dropping the last pipe-delimited field, then re-aggregate by summing frequencies.

    Parameters
    ----------
    haplotype_df : DataFrame with columns (ethnicity, haplotype, frequency)

    Returns
    -------
    DataFrame with columns (ethnicity, haplotype, frequency)
    """
    df = haplotype_df.copy()
    df['haplotype'] = df['haplotype'].str.rsplit('|', n=1).str[0]
    result = (
        df.groupby(['ethnicity', 'haplotype'], as_index=False)['frequency']
        .sum()
    )
    return result


def get_diplotype_frequencies(haplotype_df, top_k_pct=0.99):
    """
    Compute diplotype frequencies under Hardy-Weinberg equilibrium.

    Parameters
    ----------
    haplotype_df : DataFrame with columns (ethnicity, haplotype, frequency) for ONE ethnicity
    top_k_pct : float, cumulative frequency threshold; remaining haplotypes pooled as "other"

    Returns
    -------
    DataFrame with columns (haplotype1, haplotype2, diplotype_freq), sorted descending
    """
    # Sort by frequency descending
    df = haplotype_df[['haplotype', 'frequency']].copy()
    df = df.sort_values('frequency', ascending=False).reset_index(drop=True)

    # Keep top haplotypes covering top_k_pct cumulative frequency
    cumsum = df['frequency'].cumsum()
    cutoff_idx = (cumsum < top_k_pct).sum()  # number of rows before cumsum >= top_k_pct
    # Include the row that pushes us over the threshold
    cutoff_idx = min(cutoff_idx + 1, len(df))

    top = df.iloc[:cutoff_idx].copy()
    rest = df.iloc[cutoff_idx:]

    residual_freq = rest['frequency'].sum()

    haplotypes = list(top['haplotype'])
    freqs = list(top['frequency'])

    if residual_freq > 0:
        haplotypes.append('other')
        freqs.append(residual_freq)

    # Compute diplotype frequencies under HWE
    rows = []
    n = len(haplotypes)
    for i in range(n):
        for j in range(i, n):
            h1, h2 = haplotypes[i], haplotypes[j]
            f1, f2 = freqs[i], freqs[j]
            if i == j:
                diplo_freq = f1 * f1
            else:
                diplo_freq = 2.0 * f1 * f2
            rows.append({'haplotype1': h1, 'haplotype2': h2, 'diplotype_freq': diplo_freq})

    result = pd.DataFrame(rows, columns=['haplotype1', 'haplotype2', 'diplotype_freq'])
    result = result.sort_values('diplotype_freq', ascending=False).reset_index(drop=True)
    return result


def compute_coverage(diplotype_df, n_donors):
    """
    Compute expected fraction of patients who find ≥1 match in a registry of n_donors.

    Coverage(N) = Σ_g  diplotype_freq_g × [1 − (1 − diplotype_freq_g)^N]

    Parameters
    ----------
    diplotype_df : DataFrame with column 'diplotype_freq'
    n_donors : int or float, registry size

    Returns
    -------
    float in [0, 1]
    """
    freqs = np.array(diplotype_df['diplotype_freq'], dtype=np.float64)
    n = np.float64(n_donors)
    # (1 - f)^N can underflow to 0 for large N and small f — that's fine (coverage → 1 for that genotype)
    prob_no_match = np.power(1.0 - freqs, n)
    coverage = np.sum(freqs * (1.0 - prob_no_match))
    return float(np.clip(coverage, 0.0, 1.0))


def find_registry_size(diplotype_df, target_coverage, n_min=1000, n_max=10_000_000):
    """
    Binary search (log scale) for minimum N where coverage >= target_coverage.

    Parameters
    ----------
    diplotype_df : DataFrame with 'diplotype_freq' column
    target_coverage : float
    n_min : int, lower bound
    n_max : int, upper bound

    Returns
    -------
    int, minimum registry size achieving target_coverage
    """
    lo = math.log10(n_min)
    hi = math.log10(n_max)

    for _ in range(50):
        mid = (lo + hi) / 2.0
        n_mid = 10 ** mid
        cov = compute_coverage(diplotype_df, n_mid)
        if cov >= target_coverage:
            hi = mid
        else:
            lo = mid

    return math.ceil(10 ** hi)


def get_combined_haplotype_freqs(haplotype_dfs_by_ethnicity):
    """
    Compute Singapore population-weighted combined haplotype frequencies.

    Parameters
    ----------
    haplotype_dfs_by_ethnicity : dict of {ethnicity: DataFrame(haplotype, frequency)}

    Returns
    -------
    DataFrame with columns (haplotype, frequency), normalized to sum to 1.0
    """
    weighted_dfs = []
    for ethnicity, df in haplotype_dfs_by_ethnicity.items():
        weight = SG_WEIGHTS[ethnicity]
        wdf = df[['haplotype', 'frequency']].copy()
        wdf['frequency'] = wdf['frequency'] * weight
        weighted_dfs.append(wdf)

    combined = pd.concat(weighted_dfs, ignore_index=True)
    combined = combined.groupby('haplotype', as_index=False)['frequency'].sum()

    # Normalize to sum to 1.0
    total = combined['frequency'].sum()
    combined['frequency'] = combined['frequency'] / total

    return combined
