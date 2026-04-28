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

# Current BMDP+SCBB donor counts (5-locus typed) — used as vertical benchmarks
CURRENT_REGISTRY = {
    'Chinese': 45754, 'Malay': 5868, 'Indian': 5586, 'Others': 3941, 'Overall': 61149,
}


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def parse_haplotypes(haplo_df, loci, cum_freq_cap=0.995):
    """
    haplo_df: DataFrame with columns [ethnicity, haplotype, frequency]
              haplotype = pipe-separated 5-locus alleles A|B|C|DRB1|DQB1
    loci:     list of locus names (LOCI_5 or LOCI_4)
    cum_freq_cap: keep the top-K haplotypes covering this fraction of total
                  frequency mass, then renormalise to sum=1.0
    Returns:  list of (allele_tuple, frequency), sorted descending by frequency
    """
    n = len(loci)
    records = {}
    for _, row in haplo_df.iterrows():
        parts = row['haplotype'].split('|')
        key = tuple(parts[:n])
        records[key] = records.get(key, 0.0) + row['frequency']

    sorted_haps = sorted(records.items(), key=lambda x: -x[1])
    total_all = sum(v for _, v in sorted_haps)
    if total_all == 0:
        return []

    kept, cum = [], 0.0
    for k, v in sorted_haps:
        kept.append((k, v))
        cum += v / total_all
        if cum >= cum_freq_cap:
            break

    kept_total = sum(v for _, v in kept)
    return [(k, v / kept_total) for k, v in kept]


# ---------------------------------------------------------------------------
# Vectorized allele match count (integer-encoded for speed)
# ---------------------------------------------------------------------------

def allele_match_count_matrix(pat_a, pat_b, don_a, don_b):
    """
    pat_a, pat_b: shape (N_pat, n_loci) — int16 allele-index arrays
    don_a, don_b: shape (N_don, n_loci) — int16 allele-index arrays
    Returns: shape (N_pat, N_don) int8 — total allele match count per pair
    """
    pa = pat_a[:, None, :]
    pb = pat_b[:, None, :]
    da = don_a[None, :, :]
    db = don_b[None, :, :]

    assign1 = (pa == da).astype(np.int8) + (pb == db).astype(np.int8)
    assign2 = (pa == db).astype(np.int8) + (pb == da).astype(np.int8)
    locus_match = np.maximum(assign1, assign2)
    return locus_match.sum(axis=2).astype(np.int8)


def compute_partial_match_probs(haplotypes, match_thresholds, chunk_size=2000):
    """
    haplotypes:       list of (allele_tuple, frequency), normalised to sum=1
    match_thresholds: list of ints, e.g. [8, 9, 10] for 10-locus
    chunk_size:       donor diplotypes processed per iteration (caps peak RAM)
    Returns: dict {threshold: (f_diplo, p_match)} where both are np.arrays

    Alleles are integer-encoded before comparison so NumPy uses fast SIMD
    integer ops rather than slow Python-object string comparisons.
    """
    n = len(haplotypes)
    if n == 0:
        return {m: (np.array([]), np.array([])) for m in match_thresholds}

    allele_tuples = [h for h, f in haplotypes]
    freqs = np.array([f for h, f in haplotypes])

    # Encode allele strings as int16 indices for fast comparison
    all_alleles = sorted({a for h in allele_tuples for a in h})
    enc = {a: np.int16(i) for i, a in enumerate(all_alleles)}
    n_loci = len(allele_tuples[0])
    int_array = np.array(
        [[enc[a] for a in h] for h in allele_tuples], dtype=np.int16
    )  # shape (n_haps, n_loci)

    # Build diplotype index lists (upper triangle incl. diagonal)
    idx_i, idx_j, f_list = [], [], []
    for i in range(n):
        for j in range(i, n):
            idx_i.append(i)
            idx_j.append(j)
            f_list.append(freqs[i] ** 2 if i == j else 2.0 * freqs[i] * freqs[j])

    idx_i = np.array(idx_i)
    idx_j = np.array(idx_j)
    f_diplo = np.array(f_list)

    pat_a = int_array[idx_i]   # (N_diplo, n_loci)
    pat_b = int_array[idx_j]
    don_a = int_array[idx_i]
    don_b = int_array[idx_j]

    N_diplo = len(f_diplo)
    p_match_acc = {m: np.zeros(N_diplo) for m in match_thresholds}

    # Chunked donor loop — avoids allocating the full (N_diplo × N_diplo) matrix
    for start in range(0, N_diplo, chunk_size):
        end = min(start + chunk_size, N_diplo)
        chunk_match = allele_match_count_matrix(
            pat_a, pat_b, don_a[start:end], don_b[start:end]
        )  # (N_diplo, chunk_size)
        for m in match_thresholds:
            p_match_acc[m] += (chunk_match >= m).astype(np.float64) @ f_diplo[start:end]

    return {m: (f_diplo, p_match_acc[m]) for m in match_thresholds}


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

def build_combined_haplotypes(haplo_by_eth, loci, cum_freq_cap=0.95):
    """
    Weighted combination of per-ethnicity haplotypes, capped at cum_freq_cap
    of total frequency mass and renormalised to sum=1.
    Returns list of (allele_tuple, frequency), sorted descending by frequency.
    """
    combined = {}
    for eth, haps in haplo_by_eth.items():
        w = SG_WEIGHTS[eth]
        for allele_tuple, freq in haps:
            combined[allele_tuple] = combined.get(allele_tuple, 0.0) + w * freq
    total = sum(combined.values())
    if total == 0:
        return []
    sorted_haps = sorted(combined.items(), key=lambda x: -x[1])
    kept, cum = [], 0.0
    for k, v in sorted_haps:
        kept.append((k, v))
        cum += v / total
        if cum >= cum_freq_cap:
            break
    kept_total = sum(v for _, v in kept)
    return [(k, v / kept_total) for k, v in kept]


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(eth_curves_dict, match_levels, framework_name, out_path):
    """
    eth_curves_dict: {ethnicity: {min_matches: coverage_array}}
    match_levels: list of (min_matches, label, color), e.g. [(10,'10/10','blue'), ...]
                  ordered from highest to lowest
    2×2 layout: Chinese (top-left), Malay (top-right), Indian (bottom-left), Others (bottom-right)
    """
    panels = ETHNICITIES  # Chinese, Malay, Indian, Others — no Overall
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), facecolor='white')
    ax_flat = [axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]]

    for ax, eth in zip(ax_flat, panels):
        ax.set_facecolor('white')

        # Gridlines
        for y in range(10, 100, 10):
            ax.axhline(y, color='lightgrey', linestyle='--', linewidth=0.8, zorder=0)

        # Current registry size benchmark
        curr = CURRENT_REGISTRY.get(eth)
        if curr is not None:
            ax.axvline(curr, color='dimgrey', linestyle=':', linewidth=1.5,
                       zorder=1, alpha=0.85)
            ax.text(curr * 1.12, 3, f'Current\n~{curr // 1000}K',
                    fontsize=8, color='dimgrey', va='bottom', ha='left')

        curves = eth_curves_dict.get(eth, {})
        for min_m, label, color in match_levels:
            cov = curves.get(min_m, np.zeros(len(N_VALUES)))
            ax.plot(N_VALUES, cov * 100, color=color, linewidth=2.5, label=label)

        ax.set_xscale('log')
        ax.set_xlim(1, 1e6)
        ax.set_ylim(0, 100)
        ax.set_xlabel('Number of donors in the registry', fontsize=11)
        ax.set_ylabel('Percentage of patients with donors', fontsize=11)
        ax.set_title(eth, fontsize=14, fontweight='bold')
        ax.legend(fontsize=10, loc='upper left')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle(framework_name, fontsize=15, fontweight='bold', y=1.01)
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
