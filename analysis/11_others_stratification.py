"""
11_others_stratification.py
Ancestry stratification of the heterogeneous "Others" group.

Method:
  1. Pivot Others individuals to wide format (allele presence/absence across 5 loci).
  2. Build binary allele indicator matrix → PCA (n_components=10).
  3. K-means clustering on PCA scores, k=2..5, select by silhouette coefficient.
  4. Run full multi-locus EM per cluster → N* at 75/85/90/95% coverage.

Outputs:
  data/others_cluster_assignments.csv   — per-individual cluster label
  data/others_cluster_registry.csv      — N* per cluster × coverage threshold
  figures/others_pca_scatter.png        — PCA scatter coloured by cluster
  figures/others_registry_by_cluster.png — registry targets per cluster
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hwe_test import _pivot_to_individuals, _em_full_phase
from registry_model import get_diplotype_frequencies, find_registry_size

DATA_DIR = os.path.join(HERE, 'data')
FIG_DIR  = os.path.join(HERE, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

LOCI      = ['HLA-A', 'HLA-B', 'HLA-C', 'DRB1', 'DQB1']
THRESHOLDS = [0.75, 0.85, 0.90, 0.95]
FREQ_THRESHOLD = 0.001
CAP       = 5000
SEED      = 42

CLUSTER_COLORS = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00']


# ---------------------------------------------------------------------------
# Step 1-2: Build allele indicator matrix for Others
# ---------------------------------------------------------------------------

def build_allele_matrix(eth_df: pd.DataFrame, min_allele_freq: float = 0.01):
    """
    Pivot tidy df to wide, then build binary allele presence/absence matrix.
    Only includes alleles observed in >= min_allele_freq of individuals (removes outlier-creating rare alleles).
    Returns (wide_df, indicator_matrix, feature_names).
    """
    wide = _pivot_to_individuals(eth_df)
    if wide is None or len(wide) == 0:
        raise ValueError("No complete Others individuals found")

    n_indiv = len(wide)
    print(f"  Wide matrix: {n_indiv} individuals with all 5 loci typed")

    locus_cols = {loc: [f'{loc}_1', f'{loc}_2'] for loc in LOCI}

    # Count allele occurrences (each individual contributes up to 2 copies)
    allele_counts: dict = {}
    for loc, cols in locus_cols.items():
        for c in cols:
            if c not in wide.columns:
                continue
            for a in wide[c].dropna():
                allele_counts[(loc, a)] = allele_counts.get((loc, a), 0) + 1

    # Keep alleles present in >= min_allele_freq of individuals
    min_count = int(min_allele_freq * n_indiv)
    common_alleles = {(loc, a) for (loc, a), cnt in allele_counts.items() if cnt >= min_count}
    print(f"  Common alleles (≥{min_allele_freq:.0%} of individuals): {len(common_alleles)}")

    # Build indicator matrix vectorised: for each (locus, allele) column
    allele_list = sorted(common_alleles)
    col_names = [f'{loc}_{a}' for loc, a in allele_list]
    mat = np.zeros((n_indiv, len(allele_list)), dtype=np.uint8)

    for col_idx, (loc, a) in enumerate(allele_list):
        c1, c2 = f'{loc}_1', f'{loc}_2'
        present = np.zeros(n_indiv, dtype=bool)
        if c1 in wide.columns:
            present |= (wide[c1] == a).values
        if c2 in wide.columns:
            present |= (wide[c2] == a).values
        mat[:, col_idx] = present.astype(np.uint8)

    indicator = pd.DataFrame(mat, columns=col_names)
    feature_names = list(indicator.columns)
    print(f"  Allele indicator matrix: {indicator.shape[0]} × {indicator.shape[1]}")
    return wide, indicator, feature_names


# ---------------------------------------------------------------------------
# Step 3: PCA + K-means
# ---------------------------------------------------------------------------

def select_best_k(pca_scores: np.ndarray, k_range=range(2, 6)) -> tuple:
    """
    Run K-means for each k in k_range, return (best_k, labels_dict, sil_dict).
    """
    labels_dict = {}
    sil_dict    = {}
    rng = np.random.default_rng(SEED)

    for k in k_range:
        km = KMeans(n_clusters=k, n_init=20, random_state=SEED)
        labels = km.fit_predict(pca_scores)
        sil = silhouette_score(pca_scores, labels)
        labels_dict[k] = labels
        sil_dict[k] = sil
        print(f"  k={k}  silhouette={sil:.4f}")

    best_k = max(sil_dict, key=sil_dict.get)
    print(f"  → Best k={best_k} (silhouette={sil_dict[best_k]:.4f})")
    return best_k, labels_dict, sil_dict


# ---------------------------------------------------------------------------
# Step 4: EM per cluster → N*
# ---------------------------------------------------------------------------

def run_em_for_cluster(eth_df: pd.DataFrame, cluster_label: str,
                       wide: pd.DataFrame, cluster_mask: np.ndarray) -> pd.DataFrame:
    """
    Run EM on individuals in the cluster, return haplotype DataFrame.
    """
    cluster_ids = wide['sample_id'].values[cluster_mask]
    cluster_df = eth_df[eth_df['sample_id'].isin(cluster_ids)].copy()
    print(f"  Cluster {cluster_label}: {len(cluster_ids)} individuals")

    pivoted = _pivot_to_individuals(cluster_df)
    if pivoted is None or len(pivoted) == 0:
        return pd.DataFrame(columns=['haplotype', 'frequency'])

    if len(pivoted) > CAP:
        pivoted = pivoted.sample(n=CAP, random_state=SEED)

    haplo_freqs = _em_full_phase(pivoted)

    rows = [{'haplotype': '|'.join(h), 'frequency': f}
            for h, f in haplo_freqs.items() if f >= FREQ_THRESHOLD]
    df = pd.DataFrame(rows, columns=['haplotype', 'frequency'])
    df['frequency'] = df['frequency'] / df['frequency'].sum()
    print(f"    → {len(df)} haplotypes (freq ≥ {FREQ_THRESHOLD})")
    return df


def compute_registry_targets(haplo_df: pd.DataFrame) -> dict:
    """Return {threshold: N*} for a normalised haplotype DataFrame."""
    diplo = get_diplotype_frequencies(haplo_df)
    return {t: find_registry_size(diplo, t) for t in THRESHOLDS}


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_pca_scatter(pca_scores: np.ndarray, labels: np.ndarray,
                     sil_dict: dict, best_k: int, out_path: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 6), facecolor='white')
    ax.set_facecolor('white')

    ancestry_names = {0: 'Cluster 1 — European/Eurasian',
                      1: 'Cluster 2 — Filipino/SE Asian',
                      2: 'Cluster 3 — Northeast Asian/Mixed'}
    for k in range(best_k):
        mask = labels == k
        ax.scatter(pca_scores[mask, 0], pca_scores[mask, 1],
                   c=CLUSTER_COLORS[k], s=10, alpha=0.55,
                   label=ancestry_names.get(k, f'Cluster {k+1}'))
    ax.set_xlabel('PC1', fontsize=12)
    ax.set_ylabel('PC2', fontsize=12)
    ax.set_title(f'PCA of Others donors — K-means k={best_k} clusters '
                 f'(silhouette={sil_dict[best_k]:.2f})',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=10, markerscale=3, frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {out_path}')


def plot_registry_by_cluster(registry_rows: list, best_k: int, out_path: str) -> None:
    threshold_labels = {0.75: '75%', 0.85: '85%', 0.90: '90%', 0.95: '95%'}
    thresholds = THRESHOLDS
    n_clusters = best_k
    bar_width = 0.15
    x = np.arange(len(thresholds))

    fig, ax = plt.subplots(figsize=(9, 5), facecolor='white')
    ax.set_facecolor('white')

    cluster_map = {}
    for row in registry_rows:
        c = row['cluster']
        cluster_map.setdefault(c, {})
        cluster_map[c][row['target_coverage']] = row['registry_size']

    for k, (cluster, vals) in enumerate(sorted(cluster_map.items())):
        bar_vals = [vals.get(t, np.nan) for t in thresholds]
        offset = (k - (n_clusters - 1) / 2) * bar_width
        bars = ax.bar(x + offset, bar_vals, bar_width,
                      label=cluster, color=CLUSTER_COLORS[k], alpha=0.85,
                      edgecolor='white')
        for bar, v in zip(bars, bar_vals):
            if not np.isnan(v):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 200,
                        f'{int(v):,}', ha='center', va='bottom',
                        fontsize=7, rotation=90)

    ax.set_xticks(x)
    ax.set_xticklabels([threshold_labels[t] for t in thresholds], fontsize=11)
    ax.set_xlabel('Coverage target', fontsize=11)
    ax.set_ylabel('Minimum registry size (donors)', fontsize=11)
    ax.set_title('Registry size targets by Others sub-cluster\n(10/10 same-cluster matching)',
                 fontsize=11, fontweight='bold')
    ax.legend(fontsize=9)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.2)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{int(v):,}'))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {out_path}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Loading Others from hla_clean.csv...")
    hla_df = pd.read_csv(os.path.join(DATA_DIR, 'hla_clean.csv'))
    eth_df = hla_df[hla_df['ethnicity'] == 'Others'].copy()

    print("Building allele indicator matrix...")
    wide, indicator, feature_names = build_allele_matrix(eth_df)

    print("Running PCA (n_components=10)...")
    # Binary indicator matrix — no StandardScaler (rare alleles already filtered)
    pca = PCA(n_components=10, random_state=SEED)
    pca_scores = pca.fit_transform(indicator.values.astype(float))
    explained = pca.explained_variance_ratio_.cumsum()
    print(f"  Cumulative variance: PC5={explained[4]:.3f}, PC10={explained[9]:.3f}")

    print("Selecting best k by silhouette (k=2..5)...")
    best_k, labels_dict, sil_dict = select_best_k(pca_scores[:, :5])
    labels = labels_dict[best_k]

    # Persist silhouette scores so the manuscript cites a computed value,
    # not a hand-typed one (source of the v2.x silhouette=0.97 prose error).
    pd.DataFrame({'k': list(sil_dict), 'silhouette_5pc': list(sil_dict.values())}) \
        .to_csv(os.path.join(DATA_DIR, 'others_cluster_silhouette.csv'), index=False)

    # Save cluster assignments (wide has sample_id as a column after reset_index)
    assignments = pd.DataFrame({
        'sample_id': wide['sample_id'].values,
        'cluster':   [f'Cluster_{l+1}' for l in labels],
        'pc1': pca_scores[:, 0],
        'pc2': pca_scores[:, 1],
    })
    assignments.to_csv(os.path.join(DATA_DIR, 'others_cluster_assignments.csv'), index=False)
    print(f"\nCluster sizes:")
    print(assignments['cluster'].value_counts().to_string())

    print("\nRunning EM per cluster...")
    registry_rows = []
    all_haplo_rows = []

    for k in range(best_k):
        cluster_label = f'Cluster_{k+1}'
        cluster_mask = labels == k

        haplo_df = run_em_for_cluster(eth_df, cluster_label, wide, cluster_mask)
        if haplo_df.empty:
            print(f"  Skipping {cluster_label}: no haplotypes")
            continue

        haplo_df['cluster'] = cluster_label
        all_haplo_rows.append(haplo_df)

        targets = compute_registry_targets(haplo_df)
        print(f"    N* at 90%: {targets[0.90]:,}   at 95%: {targets[0.95]:,}")
        for t, n in targets.items():
            registry_rows.append({
                'cluster':         cluster_label,
                'n_individuals':   int(cluster_mask.sum()),
                'n_haplotypes':    len(haplo_df),
                'target_coverage': t,
                'registry_size':   n,
            })

    if not all_haplo_rows:
        print("ERROR: No clusters produced valid EM haplotypes. Exiting.")
        return

    # Save outputs
    haplo_all = pd.concat(all_haplo_rows, ignore_index=True)
    haplo_all.to_csv(os.path.join(DATA_DIR, 'others_cluster_haplotypes.csv'), index=False)

    registry_df = pd.DataFrame(registry_rows)
    registry_df.to_csv(os.path.join(DATA_DIR, 'others_cluster_registry.csv'), index=False)
    print(f"\nSaved: others_cluster_registry.csv")

    print("\n=== Registry size summary (10/10, 95% coverage) ===")
    mask95 = registry_df['target_coverage'] == 0.95
    print(registry_df[mask95][['cluster', 'n_individuals', 'n_haplotypes', 'registry_size']].to_string(index=False))

    print("\nGenerating figures...")
    plot_pca_scatter(pca_scores, labels, sil_dict, best_k,
                     os.path.join(FIG_DIR, 'others_pca_scatter.png'))
    plot_registry_by_cluster(registry_rows, best_k,
                              os.path.join(FIG_DIR, 'others_registry_by_cluster.png'))

    print("\nDone.")


if __name__ == '__main__':
    main()
