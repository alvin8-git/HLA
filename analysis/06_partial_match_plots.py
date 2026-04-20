"""
06_partial_match_plots.py
Partial-match registry coverage curves for CMIO populations.
Produces two figures matching WBMT Figure 5 style:
  - partial_match_10locus.png  (10/10, 9/10, 8/10)
  - partial_match_8locus.png   (8/8, 7/8, 6/8)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from itertools import combinations_with_replacement

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, 'data')
FIG_DIR = os.path.join(HERE, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# Locus definitions
LOCI_5 = ['HLA-A', 'HLA-B', 'HLA-C', 'DRB1', 'DQB1']   # 10-allele framework
LOCI_4 = ['HLA-A', 'HLA-B', 'HLA-C', 'DRB1']             # 8-allele framework

ETHNICITIES = ['Chinese', 'Malay', 'Indian', 'Others']
SG_WEIGHTS = {'Chinese': 0.77, 'Malay': 0.08, 'Indian': 0.09, 'Others': 0.06}
N_VALUES = np.logspace(3, 8, 300)


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def parse_haplotypes(haplo_df, loci):
    """
    haplo_df: DataFrame with columns [ethnicity, haplotype, frequency]
              haplotype = pipe-separated 5-locus alleles A|B|C|DRB1|DQB1
    loci:     list of locus names (LOCI_5 or LOCI_4)
    Returns:  list of (allele_tuple, frequency)
    """
    n = len(loci)  # number of loci to use (4 or 5)
    records = {}
    for _, row in haplo_df.iterrows():
        parts = row['haplotype'].split('|')
        key = tuple(parts[:n])
        records[key] = records.get(key, 0.0) + row['frequency']
    total = sum(records.values())
    if total == 0:
        return []
    return [(k, v / total) for k, v in records.items()]


def get_allele_freqs_dict(allele_freq_df, ethnicity, loci):
    """
    Returns dict: {locus: {allele: frequency}}
    Normalised per locus.
    """
    # Map locus names: allele_freqs_per_locus uses 'DRB1'/'DQB1' without prefix
    # and 'HLA-A' etc. — check actual values
    sub = allele_freq_df[allele_freq_df['ethnicity'] == ethnicity]
    result = {}
    for locus in loci:
        lsub = sub[sub['locus'] == locus]
        if lsub.empty:
            result[locus] = {}
            continue
        d = dict(zip(lsub['allele'], lsub['frequency']))
        total = sum(d.values())
        if total > 0:
            d = {a: f / total for a, f in d.items()}
        result[locus] = d
    return result


# ---------------------------------------------------------------------------
# Match probability computation
# ---------------------------------------------------------------------------

def locus_match_dist(p1, p2, freq_dict):
    """
    Returns numpy array [P(0 match), P(1 match), P(2 match)] at one locus.
    p1, p2: patient alleles at this locus
    freq_dict: {allele: frequency} for this locus
    """
    qp1 = freq_dict.get(p1, 0.0)
    qp2 = freq_dict.get(p2, 0.0)

    if p1 == p2:
        # Homozygous patient
        p = qp1
        p2m = p * p
        p1m = 2 * p * (1 - p)
        p0m = (1 - p) ** 2
    else:
        # Heterozygous patient
        p2m = 2 * qp1 * qp2
        # P(1 match): donor has exactly one of p1 or p2
        # = P(donor has p1 but not p2) + P(donor has p2 but not p1)
        # donor genotype (x,y) with x<=y under HWE:
        # More directly from multiset intersection:
        # P(1 match) = q[p1]^2 + q[p2]^2
        #            + 2*q[p1]*(1-q[p1]-q[p2])
        #            + 2*q[p2]*(1-q[p1]-q[p2])
        p1m = (qp1**2 + qp2**2
               + 2 * qp1 * (1 - qp1 - qp2)
               + 2 * qp2 * (1 - qp1 - qp2))
        p0m = 1.0 - p2m - p1m

    # Clip to [0, 1] to handle floating point edge cases
    arr = np.array([p0m, p1m, p2m], dtype=np.float64)
    arr = np.clip(arr, 0.0, 1.0)
    arr /= arr.sum()
    return arr


def compute_diplotype_partial_match_probs(haplotypes, allele_freqs, loci, min_matches):
    """
    For each diplotype (h_p, h_q):
      - compute P(total allele matches >= min_matches) against a random donor
    Returns: (f_g_arr, p_match_arr) as numpy arrays.
    """
    n_loci = len(loci)
    max_matches = 2 * n_loci

    f_g_list = []
    p_match_list = []

    haps = haplotypes  # list of (allele_tuple, freq)
    n_haps = len(haps)

    for i in range(n_haps):
        h_p, f_p = haps[i]
        for j in range(i, n_haps):
            h_q, f_q = haps[j]

            # Diplotype frequency under HWE
            if i == j:
                f_g = f_p * f_q          # homozygous diplotype
            else:
                f_g = 2 * f_p * f_q      # heterozygous

            if f_g < 1e-12:
                continue

            # Convolve locus match distributions
            d = np.zeros(max_matches + 1)
            d[0] = 1.0

            for li, locus in enumerate(loci):
                p1 = h_p[li]
                p2 = h_q[li]
                lmd = locus_match_dist(p1, p2, allele_freqs[locus])
                d = np.convolve(d, lmd)[:max_matches + 1]

            # P(total >= min_matches)
            p_m = float(d[min_matches:].sum())

            f_g_list.append(f_g)
            p_match_list.append(p_m)

    return np.array(f_g_list), np.array(p_match_list)


def coverage_curve(f_g_arr, p_m_arr, n_values):
    """
    Vectorised coverage curve computation.
    Returns array of shape (len(n_values),).
    """
    if len(f_g_arr) == 0:
        return np.zeros(len(n_values))
    # Shape: (K, len(N)) — then sum over K
    contrib = f_g_arr[:, None] * (1 - np.power(
        np.clip(1 - p_m_arr[:, None], 0, 1),
        n_values[None, :]
    ))
    return contrib.sum(axis=0)


# ---------------------------------------------------------------------------
# Combined (Overall) haplotype pool
# ---------------------------------------------------------------------------

def build_combined_haplotypes(haplo_by_eth, loci):
    """
    Weighted combination of per-ethnicity haplotypes.
    Returns list of (allele_tuple, frequency).
    """
    combined = {}
    for eth, haps in haplo_by_eth.items():
        w = SG_WEIGHTS[eth]
        for allele_tuple, freq in haps:
            combined[allele_tuple] = combined.get(allele_tuple, 0.0) + w * freq
    total = sum(combined.values())
    if total == 0:
        return []
    return [(k, v / total) for k, v in combined.items()]


def build_combined_allele_freqs(allele_freqs_by_eth, loci):
    """
    Weighted average allele frequencies across ethnicities.
    Returns {locus: {allele: freq}}, normalised.
    """
    combined = {locus: {} for locus in loci}
    for eth, afreqs in allele_freqs_by_eth.items():
        w = SG_WEIGHTS[eth]
        for locus in loci:
            for allele, freq in afreqs.get(locus, {}).items():
                combined[locus][allele] = combined[locus].get(allele, 0.0) + w * freq
    # Normalise each locus
    for locus in loci:
        total = sum(combined[locus].values())
        if total > 0:
            combined[locus] = {a: f / total for a, f in combined[locus].items()}
    return combined


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(eth_curves_dict, match_levels, framework_name, out_path):
    """
    eth_curves_dict: {ethnicity: {min_matches: coverage_array}}
    match_levels: list of (min_matches, label, color), e.g. [(10,'10/10','blue'), ...]
                  ordered from highest to lowest
    """
    panels = ETHNICITIES + ['Overall']
    fig, axes = plt.subplots(1, 5, figsize=(25, 5), facecolor='white')

    for ax, eth in zip(axes, panels):
        ax.set_facecolor('white')

        # Gridlines
        for y in range(10, 100, 10):
            ax.axhline(y, color='lightgrey', linestyle='--', linewidth=0.8, zorder=0)

        curves = eth_curves_dict.get(eth, {})
        for min_m, label, color in match_levels:
            cov = curves.get(min_m, np.zeros(len(N_VALUES)))
            ax.plot(N_VALUES, cov * 100, color=color, linewidth=2, label=label)

        ax.set_xscale('log')
        ax.set_xlim(1e3, 1e8)
        ax.set_ylim(0, 100)
        ax.set_xlabel('Number of donors in the registry', fontsize=10)
        ax.set_ylabel('Percentage of patients with donors', fontsize=10)
        ax.set_title(eth, fontsize=12, fontweight='bold')
        ax.legend(fontsize=9, loc='upper left')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle(framework_name, fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # Load data
    haplo_df = pd.read_csv(os.path.join(DATA_DIR, 'haplo_freqs_em.csv'))
    allele_df = pd.read_csv(os.path.join(DATA_DIR, 'allele_freqs_per_locus.csv'))

    # Check locus names in allele_df
    print("Loci in allele_freqs:", sorted(allele_df['locus'].unique()))

    for framework_name, loci, match_levels, out_fname in [
        (
            '10-locus (HLA-A, B, C, DRB1, DQB1) matching',
            LOCI_5,
            [(10, '10/10', 'blue'), (9, '9/10', 'red'), (8, '8/10', 'green')],
            'partial_match_10locus.png',
        ),
        (
            '8-locus (HLA-A, B, C, DRB1) matching',
            LOCI_4,
            [(8, '8/8', 'blue'), (7, '7/8', 'red'), (6, '6/8', 'green')],
            'partial_match_8locus.png',
        ),
    ]:
        print(f"\n=== {framework_name} ===")
        eth_curves_dict = {}
        haplo_by_eth = {}
        afreqs_by_eth = {}

        for eth in ETHNICITIES:
            print(f"  Processing {eth}...")
            eth_haplo_df = haplo_df[haplo_df['ethnicity'] == eth]
            haplotypes = parse_haplotypes(eth_haplo_df, loci)
            allele_freqs = get_allele_freqs_dict(allele_df, eth, loci)
            haplo_by_eth[eth] = haplotypes
            afreqs_by_eth[eth] = allele_freqs

            curves = {}
            for min_m, label, color in match_levels:
                f_g_arr, p_m_arr = compute_diplotype_partial_match_probs(
                    haplotypes, allele_freqs, loci, min_m
                )
                curves[min_m] = coverage_curve(f_g_arr, p_m_arr, N_VALUES)
            eth_curves_dict[eth] = curves

            # Sanity check at N=1,000,000
            idx_1m = np.searchsorted(N_VALUES, 1e6)
            for min_m, label, _ in match_levels:
                val = curves[min_m][idx_1m] * 100
                print(f"    {label} @ N=1M: {val:.1f}%")

        # Overall
        print("  Processing Overall...")
        combined_haps = build_combined_haplotypes(haplo_by_eth, loci)
        combined_afreqs = build_combined_allele_freqs(afreqs_by_eth, loci)
        curves = {}
        for min_m, label, color in match_levels:
            f_g_arr, p_m_arr = compute_diplotype_partial_match_probs(
                combined_haps, combined_afreqs, loci, min_m
            )
            curves[min_m] = coverage_curve(f_g_arr, p_m_arr, N_VALUES)
        eth_curves_dict['Overall'] = curves

        idx_1m = np.searchsorted(N_VALUES, 1e6)
        for min_m, label, _ in match_levels:
            val = curves[min_m][idx_1m] * 100
            print(f"    Overall {label} @ N=1M: {val:.1f}%")

        out_path = os.path.join(FIG_DIR, out_fname)
        make_figure(eth_curves_dict, match_levels, framework_name, out_path)

    print("\nDone.")


if __name__ == '__main__':
    main()
