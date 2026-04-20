"""
06_partial_match_plots.py
Partial-match registry coverage curves for CMIO populations.
Produces two figures matching WBMT Figure 5 style:
  - partial_match_10locus.png  (10/10, 9/10, 8/10)
  - partial_match_8locus.png   (8/8, 7/8, 6/8)

Uses haplotype-pair enumeration to correctly capture linkage disequilibrium (LD).
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, 'data')
FIG_DIR = os.path.join(HERE, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# Locus definitions
LOCI_5 = ['HLA-A', 'HLA-B', 'HLA-C', 'DRB1', 'DQB1']   # 10-allele framework
LOCI_4 = ['HLA-A', 'HLA-B', 'HLA-C', 'DRB1']             # 8-allele framework

ETHNICITIES = ['Chinese', 'Malay', 'Indian', 'Others']
SG_WEIGHTS = {'Chinese': 0.77, 'Malay': 0.08, 'Indian': 0.09, 'Others': 0.06}
N_VALUES = np.logspace(0, 6, 300)


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def parse_haplotypes(haplo_df, loci):
    """
    haplo_df: DataFrame with columns [ethnicity, haplotype, frequency]
              haplotype = pipe-separated 5-locus alleles A|B|C|DRB1|DQB1
    loci:     list of locus names (LOCI_5 or LOCI_4)
    Returns:  list of (allele_tuple, frequency), normalized to sum=1.0
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
    # Normalize so frequencies sum to 1.0 (EM output may not include all haplotypes)
    return [(k, v / total) for k, v in records.items()]


# ---------------------------------------------------------------------------
# Vectorized allele match count
# ---------------------------------------------------------------------------

def allele_match_count_matrix(pat_a, pat_b, don_a, don_b):
    """
    pat_a, pat_b: shape (N_pat, n_loci) — patient allele arrays (object dtype)
    don_a, don_b: shape (N_don, n_loci) — donor allele arrays (object dtype)
    Returns: shape (N_pat, N_don) int8 — total allele match count per pair
    """
    # Expand dims for broadcasting: (N_pat, 1, n_loci) vs (1, N_don, n_loci)
    pa = pat_a[:, None, :]   # (N_pat, 1, n_loci)
    pb = pat_b[:, None, :]
    da = don_a[None, :, :]   # (1, N_don, n_loci)
    db = don_b[None, :, :]

    # Two possible allele assignments at each locus, take the better one
    assign1 = (pa == da).astype(np.int8) + (pb == db).astype(np.int8)
    assign2 = (pa == db).astype(np.int8) + (pb == da).astype(np.int8)
    locus_match = np.maximum(assign1, assign2)  # (N_pat, N_don, n_loci)
    return locus_match.sum(axis=2).astype(np.int8)  # (N_pat, N_don)


def compute_partial_match_probs(haplotypes, match_thresholds):
    """
    haplotypes: list of (allele_tuple, frequency), normalized to sum=1
    match_thresholds: list of ints, e.g. [8, 9, 10] for 10-locus
    Returns: dict {threshold: (f_diplo, p_match)} where both are np.arrays
    """
    n = len(haplotypes)
    if n == 0:
        return {m: (np.array([]), np.array([])) for m in match_thresholds}

    allele_tuples = [h for h, f in haplotypes]
    freqs = np.array([f for h, f in haplotypes])

    # Parse haplotype alleles into arrays: shape (n, n_loci)
    allele_array = np.array(allele_tuples, dtype=object)  # (n, n_loci)

    # Build diplotype list (upper triangle including diagonal)
    diplo_list = []
    for i in range(n):
        for j in range(i, n):
            f = freqs[i] ** 2 if i == j else 2 * freqs[i] * freqs[j]
            diplo_list.append((i, j, f))

    idx_i = np.array([d[0] for d in diplo_list])
    idx_j = np.array([d[1] for d in diplo_list])
    f_diplo = np.array([d[2] for d in diplo_list])

    # Patient and donor diplotype allele arrays (same population pool)
    pat_a = allele_array[idx_i]   # (N_diplo, n_loci)
    pat_b = allele_array[idx_j]
    don_a = allele_array[idx_i]   # (N_diplo, n_loci)
    don_b = allele_array[idx_j]

    # Compute allele match count matrix: (N_diplo_patient, N_diplo_donor)
    total_match = allele_match_count_matrix(pat_a, pat_b, don_a, don_b)

    results = {}
    for m in match_thresholds:
        # p_match[i] = Σ_j f_diplo[j] * I[total_match[i,j] >= m]
        p_match = (total_match >= m).astype(np.float64) @ f_diplo
        results[m] = (f_diplo, p_match)

    return results


def coverage_curve(f_diplo, p_match, n_values):
    """
    Vectorised coverage curve computation.
    Returns array of shape (len(n_values),).
    """
    if len(f_diplo) == 0:
        return np.zeros(len(n_values))
    # Shape: (K, len(N)) — then sum over K
    contrib = f_diplo[:, None] * (1.0 - np.power(
        np.clip(1.0 - p_match[:, None], 0, 1),
        n_values[None, :]
    ))
    return contrib.sum(axis=0)


# ---------------------------------------------------------------------------
# Combined (Overall) haplotype pool
# ---------------------------------------------------------------------------

def build_combined_haplotypes(haplo_by_eth, loci):
    """
    Weighted combination of per-ethnicity haplotypes.
    Returns list of (allele_tuple, frequency), normalized to sum=1.
    """
    combined = {}
    for eth, haps in haplo_by_eth.items():
        w = SG_WEIGHTS[eth]
        for allele_tuple, freq in haps:
            combined[allele_tuple] = combined.get(allele_tuple, 0.0) + w * freq
    total = sum(combined.values())
    if total == 0:
        return []
    return [(k, v / total) for k, v in sorted(combined.items(), key=lambda x: -x[1])]


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
        ax.set_xlim(1, 1e6)
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

    for framework_name, loci, match_levels, out_fname in [
        (
            '10-locus (HLA-A, B, C, DRB1, DQB1) matching',
            LOCI_5,
            [(10, '10/10', 'green'), (9, '9/10', 'blue'), (8, '8/10', 'red')],
            'partial_match_10locus.png',
        ),
        (
            '8-locus (HLA-A, B, C, DRB1) matching',
            LOCI_4,
            [(8, '8/8', 'green'), (7, '7/8', 'blue'), (6, '6/8', 'red')],
            'partial_match_8locus.png',
        ),
    ]:
        print(f"\n=== {framework_name} ===")
        eth_curves_dict = {}
        haplo_by_eth = {}
        match_thresholds = [ml[0] for ml in match_levels]

        for eth in ETHNICITIES:
            print(f"  Processing {eth}...")
            eth_haplo_df = haplo_df[haplo_df['ethnicity'] == eth]
            haplotypes = parse_haplotypes(eth_haplo_df, loci)
            haplo_by_eth[eth] = haplotypes
            print(f"    {len(haplotypes)} haplotypes loaded")

            partial_match_results = compute_partial_match_probs(haplotypes, match_thresholds)

            curves = {}
            for min_m, label, color in match_levels:
                f_diplo, p_match = partial_match_results[min_m]
                curves[min_m] = coverage_curve(f_diplo, p_match, N_VALUES)
            eth_curves_dict[eth] = curves

            # Sanity check at N=100K and N=1M
            idx_100k = np.searchsorted(N_VALUES, 1e5)
            idx_1m = np.searchsorted(N_VALUES, 1e6)
            for min_m, label, _ in match_levels:
                v100k = curves[min_m][idx_100k] * 100
                v1m = curves[min_m][idx_1m] * 100
                print(f"    {label} @ N=100K: {v100k:.1f}%   @ N=1M: {v1m:.1f}%")

        # Overall
        print("  Processing Overall...")
        combined_haps = build_combined_haplotypes(haplo_by_eth, loci)
        print(f"    {len(combined_haps)} combined haplotypes")
        partial_match_results = compute_partial_match_probs(combined_haps, match_thresholds)
        curves = {}
        for min_m, label, color in match_levels:
            f_diplo, p_match = partial_match_results[min_m]
            curves[min_m] = coverage_curve(f_diplo, p_match, N_VALUES)
        eth_curves_dict['Overall'] = curves

        idx_100k = np.searchsorted(N_VALUES, 1e5)
        idx_1m = np.searchsorted(N_VALUES, 1e6)
        for min_m, label, _ in match_levels:
            v100k = curves[min_m][idx_100k] * 100
            v1m = curves[min_m][idx_1m] * 100
            print(f"    Overall {label} @ N=100K: {v100k:.1f}%   @ N=1M: {v1m:.1f}%")

        out_path = os.path.join(FIG_DIR, out_fname)
        make_figure(eth_curves_dict, match_levels, framework_name, out_path)

    print("\nDone.")


if __name__ == '__main__':
    main()
