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


def run_em_haplotypes(df: pd.DataFrame, max_iter: int = 200,
                      tol: float = 1e-6, freq_threshold: float = 0.000001,
                      # cap raised 5,000 -> 50,000 (2026-08-20). At the 1e-3 floor
                      # the cap cost ~8%; at 1e-4 it inflates Chinese N* by 264%
                      # (11,487,962 capped vs 3,153,571 at the full 45,018 sample),
                      # because a 5,000-individual EM cannot resolve the rare tail
                      # and retains spurious phase-ambiguity haplotypes. 50,000
                      # exceeds every CMIO group, so the cap no longer binds.
                      cap: int = 50000, random_state: int = 42) -> pd.DataFrame:
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

        haplo_freqs = _em_full_phase(pivoted, max_iter=max_iter, tol=tol)

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


def _enum_phase_configs(geno: list) -> list:
    """
    Enumerate all distinct phase configurations for a multi-locus genotype.

    Parameters
    ----------
    geno : list of (a1, a2) pairs, one per locus.
           a1 == a2 indicates a homozygous locus.

    Returns
    -------
    List of (h1_tuple, h2_tuple) pairs.  The first heterozygous locus is
    always assigned a1→h1 (canonical form), so each unordered diplotype
    appears exactly once.  A fully homozygous genotype returns one entry.
    """
    het_idx = [i for i, (a1, a2) in enumerate(geno) if a1 != a2]
    n_het = len(het_idx)

    base_h1 = [a1 for a1, _ in geno]
    base_h2 = [a2 for _, a2 in geno]

    if n_het <= 1:
        return [(tuple(base_h1), tuple(base_h2))]

    configs = []
    # Fix first het locus on h1 (removes h1↔h2 symmetry).
    # Vary the remaining n_het-1 loci over all 2^(n_het-1) assignments.
    for bits in range(2 ** (n_het - 1)):
        h1 = base_h1[:]
        h2 = base_h2[:]
        for j, locus_i in enumerate(het_idx[1:]):
            if (bits >> j) & 1:
                h1[locus_i], h2[locus_i] = h2[locus_i], h1[locus_i]
        configs.append((tuple(h1), tuple(h2)))
    return configs


def _em_full_phase(wide: pd.DataFrame, max_iter: int = 200,
                   tol: float = 1e-6) -> dict:
    """
    Full multi-locus EM for 5-locus haplotype frequency estimation.

    Unlike the previous product-approximation approach, this correctly
    handles phase uncertainty: for each individual the E-step weights every
    valid phase configuration by f(h1)*f(h2) under the current frequency
    estimates, then the M-step accumulates fractional haplotype counts.

    Returns dict: haplotype_tuple → frequency
    """
    locus_cols1 = [f"{loc}_1" for loc in LOCI]
    locus_cols2 = [f"{loc}_2" for loc in LOCI]

    # Build genotype as list of (a1, a2) per locus; treat missing allele2
    # as homozygous (allele2 = allele1).
    genotypes = []
    for _, row in wide.iterrows():
        geno = []
        for c1, c2 in zip(locus_cols1, locus_cols2):
            a1 = str(row[c1]) if pd.notna(row[c1]) else None
            a2 = str(row[c2]) if pd.notna(row[c2]) else a1
            if a1 is None:
                break
            geno.append((a1, a2))
        if len(geno) == len(LOCI):
            genotypes.append(geno)

    if not genotypes:
        return {}

    # Pre-enumerate phase configs per individual (done once, outside EM loop)
    individual_configs = [_enum_phase_configs(g) for g in genotypes]

    # Collect all haplotypes that appear in any valid config
    all_haps = list({h for configs in individual_configs
                     for h1, h2 in configs for h in (h1, h2)})
    hap_idx = {h: i for i, h in enumerate(all_haps)}
    n_haps = len(all_haps)

    freqs = np.ones(n_haps) / n_haps  # uniform initialisation

    for _ in range(max_iter):
        counts = np.zeros(n_haps)

        for configs in individual_configs:
            if len(configs) == 1:
                # No phase ambiguity — contribute directly
                h1, h2 = configs[0]
                i1, i2 = hap_idx[h1], hap_idx[h2]
                if i1 == i2:
                    counts[i1] += 2.0
                else:
                    counts[i1] += 1.0
                    counts[i2] += 1.0
            else:
                # E-step: weight each config by diplotype frequency
                weights = np.empty(len(configs))
                for k, (h1, h2) in enumerate(configs):
                    i1, i2 = hap_idx[h1], hap_idx[h2]
                    weights[k] = (freqs[i1] ** 2 if i1 == i2
                                  else 2.0 * freqs[i1] * freqs[i2])

                total_w = weights.sum()
                if total_w == 0.0:
                    weights = np.ones(len(configs)) / len(configs)
                    total_w = 1.0

                # M-step: accumulate fractional counts
                for (h1, h2), w in zip(configs, weights):
                    wn = w / total_w
                    i1, i2 = hap_idx[h1], hap_idx[h2]
                    if i1 == i2:
                        counts[i1] += 2.0 * wn
                    else:
                        counts[i1] += wn
                        counts[i2] += wn

        total = counts.sum()
        new_freqs = counts / total if total > 0 else freqs
        delta = float(np.max(np.abs(new_freqs - freqs)))
        freqs = new_freqs
        if delta < tol:
            break

    return {all_haps[i]: float(freqs[i]) for i in range(n_haps)}
