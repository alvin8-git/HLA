"""
10_ld_report.py
Pairwise linkage disequilibrium (D', r²) between the 5 HLA loci,
computed from EM haplotype frequencies for each CMIO ethnicity.

For each locus pair (L1, L2) and each ethnicity, we compute:
  - Composite D' = Σ_{a,b} p_a q_b |D'_{ab}|  (frequency-weighted mean |D'|)
  - Composite r²  = Σ_{a,b} p_a q_b  r²_{ab}  (frequency-weighted mean r²)

These composite measures summarise multi-allelic LD as a single number per
locus pair, as recommended by Garner & Slatkin (2003).

Outputs:
  data/ld_report.csv           — composite D' and r² per (ethnicity, locus pair)
  figures/ld_heatmap_dprime.png — heatmap of composite D' (all ethnicities)
  figures/ld_heatmap_r2.png    — heatmap of composite r² (all ethnicities)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

HERE    = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, 'data')
FIG_DIR  = os.path.join(HERE, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

LOCI_5     = ['HLA-A', 'HLA-B', 'HLA-C', 'DRB1', 'DQB1']
LOCI_SHORT = ['A', 'B', 'C', 'DRB1', 'DQB1']
ETHNICITIES = ['Chinese', 'Malay', 'Indian', 'Others']


# ---------------------------------------------------------------------------
# Core LD computation
# ---------------------------------------------------------------------------

def compute_ld_pair(haplotypes: list, freqs: np.ndarray,
                    idx1: int, idx2: int) -> tuple:
    """
    Compute composite D' and r² for one locus pair (idx1, idx2).

    haplotypes : list of tuples, one tuple per haplotype (5-locus alleles)
    freqs      : 1-D array of haplotype frequencies (normalised, sum=1)
    idx1, idx2 : locus indices within the 5-locus tuple

    Returns (composite_dprime, composite_r2)
    """
    # Build joint allele pair frequency table f(a, b) = sum of haplotype freqs
    # where haplotype has allele a at locus idx1 and allele b at locus idx2
    joint: dict = {}
    marg1: dict = {}
    marg2: dict = {}

    for hap, f in zip(haplotypes, freqs):
        a = hap[idx1]
        b = hap[idx2]
        joint[(a, b)]  = joint.get((a, b), 0.0) + f
        marg1[a]       = marg1.get(a, 0.0) + f
        marg2[b]       = marg2.get(b, 0.0) + f

    alleles1 = list(marg1.keys())
    alleles2 = list(marg2.keys())

    comp_dprime = 0.0
    comp_r2     = 0.0

    for a in alleles1:
        p_a = marg1[a]
        for b in alleles2:
            p_b   = marg2[b]
            f_ab  = joint.get((a, b), 0.0)
            D     = f_ab - p_a * p_b

            # D_max
            if D >= 0:
                D_max = min(p_a * (1 - p_b), (1 - p_a) * p_b)
            else:
                D_max = max(-p_a * p_b, -(1 - p_a) * (1 - p_b))

            if abs(D_max) < 1e-10:
                continue

            dprime = D / D_max
            denom_r2 = p_a * (1 - p_a) * p_b * (1 - p_b)
            r2 = D ** 2 / denom_r2 if denom_r2 > 1e-10 else 0.0

            weight = p_a * p_b
            comp_dprime += weight * abs(dprime)
            comp_r2     += weight * r2

    return float(comp_dprime), float(comp_r2)


def compute_ld_all_pairs(haplo_df: pd.DataFrame, eth: str) -> list:
    """Return list of dicts with LD stats for all 10 locus pairs."""
    sub = haplo_df[haplo_df['ethnicity'] == eth].copy()
    freqs = sub['frequency'].values
    freqs = freqs / freqs.sum()

    # Parse haplotype strings to tuples
    haplotypes = [tuple(h.split('|')) for h in sub['haplotype']]

    rows = []
    for i in range(len(LOCI_5)):
        for j in range(i + 1, len(LOCI_5)):
            dprime, r2 = compute_ld_pair(haplotypes, freqs, i, j)
            rows.append({
                'ethnicity': eth,
                'locus1': LOCI_5[i],
                'locus2': LOCI_5[j],
                'locus_pair': f'{LOCI_SHORT[i]}–{LOCI_SHORT[j]}',
                'composite_dprime': round(dprime, 4),
                'composite_r2': round(r2, 4),
            })
    return rows


# ---------------------------------------------------------------------------
# Heatmap figures
# ---------------------------------------------------------------------------

def make_ld_heatmap(ld_df: pd.DataFrame, metric: str, out_path: str,
                    title: str, cmap: str) -> None:
    """4-panel heatmap (one per ethnicity)."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), facecolor='white')
    axes = axes.flatten()

    for ax, eth in zip(axes, ETHNICITIES):
        sub = ld_df[ld_df['ethnicity'] == eth]
        # Build symmetric 5×5 matrix
        mat = pd.DataFrame(np.nan, index=LOCI_SHORT, columns=LOCI_SHORT)
        for _, row in sub.iterrows():
            l1 = LOCI_SHORT[LOCI_5.index(row['locus1'])]
            l2 = LOCI_SHORT[LOCI_5.index(row['locus2'])]
            mat.loc[l1, l2] = row[metric]
            mat.loc[l2, l1] = row[metric]
        for locus in LOCI_SHORT:
            mat.loc[locus, locus] = 1.0

        sns.heatmap(mat.astype(float), ax=ax, cmap=cmap, vmin=0, vmax=1,
                    annot=True, fmt='.2f', linewidths=0.5,
                    cbar_kws={'shrink': 0.8}, square=True)
        ax.set_title(eth, fontsize=12, fontweight='bold')
        ax.tick_params(axis='both', labelsize=10)

    fig.suptitle(title, fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {out_path}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    haplo_df = pd.read_csv(os.path.join(DATA_DIR, 'haplo_freqs_em.csv'))

    all_rows = []
    for eth in ETHNICITIES:
        print(f'  Computing LD for {eth}...')
        rows = compute_ld_all_pairs(haplo_df, eth)
        all_rows.extend(rows)
        for r in rows:
            print(f"    {r['locus_pair']:12s}  D'={r['composite_dprime']:.3f}  r²={r['composite_r2']:.3f}")

    ld_df = pd.DataFrame(all_rows)
    out_csv = os.path.join(DATA_DIR, 'ld_report.csv')
    ld_df.to_csv(out_csv, index=False)
    print(f'\nSaved: {out_csv}')

    print('\nGenerating D\' heatmap...')
    make_ld_heatmap(
        ld_df, 'composite_dprime',
        os.path.join(FIG_DIR, 'ld_heatmap_dprime.png'),
        "Composite D' between HLA loci (frequency-weighted mean |D'|)",
        'YlOrRd',
    )

    print('Generating r² heatmap...')
    make_ld_heatmap(
        ld_df, 'composite_r2',
        os.path.join(FIG_DIR, 'ld_heatmap_r2.png'),
        'Composite r² between HLA loci (frequency-weighted mean r²)',
        'Blues',
    )

    print('\nDone.')


if __name__ == '__main__':
    main()
