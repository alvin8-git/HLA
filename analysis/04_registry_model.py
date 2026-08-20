"""
04_registry_model.py — Full registry model analysis (IMPL-6).
Generates coverage curves and registry size targets for HLA matching.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from registry_model import (
    collapse_to_4locus,
    get_diplotype_frequencies,
    compute_coverage,
    find_registry_size,
    get_combined_haplotype_freqs,
    diplotype_freq_vector,
    find_registry_size_vec,
    SG_WEIGHTS,
)

ETHNICITIES = ['Chinese', 'Malay', 'Indian', 'Others']
THRESHOLDS = [0.75, 0.85, 0.90, 0.95]
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, 'data')
FIG_DIR = os.path.join(HERE, 'figures')  # was HERE/../figures; nothing reads that dir

N_SWEEP = np.logspace(3, 9, 200)


_CROSS_ALIGN_CACHE = {}


def _donor_kept_map(donor_haplo_df, top_k_pct=0.99):
    """
    Reproduce get_diplotype_frequencies' kept-set + 'other' pooling as a
    {haplotype: frequency} lookup, WITHOUT expanding to diplotypes.
    """
    d = donor_haplo_df[['haplotype', 'frequency']].copy()
    d = d.sort_values('frequency', ascending=False).reset_index(drop=True)
    d['frequency'] = d['frequency'] / d['frequency'].sum()
    cut = min(int((d['frequency'].cumsum() < top_k_pct).sum()) + 1, len(d))
    m = dict(zip(d['haplotype'][:cut], d['frequency'][:cut]))
    residual = float(d['frequency'][cut:].sum())
    if residual > 0:
        m['other'] = residual
    return m


def _align_cross(patient_diplo, donor_haplo_df):
    """
    Align patient and donor diplotype frequency vectors once, and cache.

    The donor side is a haplotype-level lookup rather than a diplotype merge.
    At a 1e-6 floor the combined donor pool is ~25,000 haplotypes, i.e. ~300M
    diplotype pairs; materialising that frame (two object columns) to merge
    against would exhaust memory, and the merge result never depends on N
    anyway. Because get_diplotype_frequencies now labels pairs in canonical
    lexicographic order, looking the two haplotypes up independently and
    recombining under HWE reproduces the merge exactly (verified to 1.7e-18).
    """
    key = (id(patient_diplo), id(donor_haplo_df))
    hit = _CROSS_ALIGN_CACHE.get(key)
    if hit is not None:
        return hit
    cmap = _donor_kept_map(donor_haplo_df)
    h1 = patient_diplo['haplotype1'].to_numpy()
    h2 = patient_diplo['haplotype2'].to_numpy()
    g1 = np.fromiter((cmap.get(x, 0.0) for x in h1), dtype=np.float64, count=len(h1))
    g2 = np.fromiter((cmap.get(x, 0.0) for x in h2), dtype=np.float64, count=len(h2))
    donor_freq = np.where(h1 == h2, g1 * g2, 2.0 * g1 * g2)
    pf = patient_diplo['diplotype_freq'].to_numpy(dtype=np.float64)
    _CROSS_ALIGN_CACHE.clear()          # only ever need the current pair
    _CROSS_ALIGN_CACHE[key] = (pf, donor_freq)
    return pf, donor_freq


def compute_coverage_cross(patient_diplo, donor_haplo_df, n_donors):
    """
    Coverage where patients are drawn from patient_diplo and donors from the
    haplotype pool donor_haplo_df.

    Coverage_cross(N) = Σ_g  fg_patient × [1 − (1 − fg_donor)^N]
    """
    pf, donor_freq = _align_cross(patient_diplo, donor_haplo_df)
    p_match = 1.0 - np.power(1.0 - donor_freq, float(n_donors))
    return float(np.clip((pf * p_match).sum(), 0.0, 1.0))


def find_registry_size_cross(patient_diplo, donor_diplo, target_coverage,
                              n_min=1000, n_max=10_000_000_000):
    """Binary search on log scale for cross-ethnic coverage."""
    import math
    lo = math.log10(n_min)
    hi = math.log10(n_max)
    for _ in range(50):
        mid = (lo + hi) / 2.0
        n_mid = 10 ** mid
        cov = compute_coverage_cross(patient_diplo, donor_diplo, n_mid)
        if cov >= target_coverage:
            hi = mid
        else:
            lo = mid
    return math.ceil(10 ** hi)


def run_sweep(patient_diplo, donor_diplo, n_sweep, cross=False):
    """Return array of coverage values for each N in n_sweep."""
    coverages = []
    for n in n_sweep:
        if cross:
            cov = compute_coverage_cross(patient_diplo, donor_diplo, n)
        else:
            cov = compute_coverage(patient_diplo, n)
        coverages.append(cov)
    return np.array(coverages)


def make_figure(curve_rows, match_level, fig_path):
    """
    Plot coverage curves for all ethnicities including Combined.

    curve_rows: list of dicts with keys:
        ethnicity, model_variant, N (array), coverage (array)
    """
    all_ethnicities = ETHNICITIES + ['Combined']
    fig, axes = plt.subplots(1, 5, figsize=(20, 4), sharey=True)

    colors = {'same_ethnicity': 'steelblue', 'cross_ethnic': 'darkorange'}
    styles = {'same_ethnicity': '-', 'cross_ethnic': '--'}
    labels = {'same_ethnicity': 'Same ethnicity', 'cross_ethnic': 'Cross-ethnic (SG pool)'}

    for ax, eth in zip(axes, all_ethnicities):
        eth_rows = [r for r in curve_rows if r['ethnicity'] == eth]
        for row in eth_rows:
            mv = row['model_variant']
            ax.plot(row['N'], row['coverage'],
                    color=colors[mv], linestyle=styles[mv],
                    label=labels[mv], linewidth=1.5)

        for threshold in THRESHOLDS:
            ax.axhline(threshold, color='grey', linestyle='--', linewidth=0.8, alpha=0.6)

        ax.set_xscale('log')
        ax.set_xlim(N_SWEEP[0], N_SWEEP[-1])
        ax.set_ylim(0, 1)
        ax.set_title(eth)
        ax.set_xlabel('Registry size (N)')

    axes[0].set_ylabel('Coverage')

    # Legend in first panel only
    handles, lbls = axes[0].get_legend_handles_labels()
    # Deduplicate
    seen = {}
    for h, l in zip(handles, lbls):
        if l not in seen:
            seen[l] = h
    axes[0].legend(seen.values(), seen.keys(), fontsize=8, loc='lower right')

    fig.suptitle(f'Registry Coverage — {match_level} match', fontsize=13, y=1.01)
    plt.tight_layout()
    os.makedirs(os.path.dirname(fig_path), exist_ok=True)
    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved figure: {fig_path}")


def main():
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)

    # 1. Load haplotype frequencies
    haplo5 = pd.read_csv(os.path.join(DATA_DIR, 'haplo_freqs_em.csv'))
    haplo4 = collapse_to_4locus(haplo5)

    # 2. Build per-ethnicity haplotype dfs (normalize within each ethnicity so freqs sum to 1)
    def extract_and_normalize(df, eth):
        sub = df[df['ethnicity'] == eth][['haplotype', 'frequency']].copy().reset_index(drop=True)
        sub['frequency'] = sub['frequency'] / sub['frequency'].sum()
        return sub

    haplo5_by_eth = {eth: extract_and_normalize(haplo5, eth) for eth in ETHNICITIES}
    haplo4_by_eth = {eth: extract_and_normalize(haplo4, eth) for eth in ETHNICITIES}

    # 3. Combined haplotype freqs (SG-weighted)
    combined5 = get_combined_haplotype_freqs(haplo5_by_eth)
    combined4 = get_combined_haplotype_freqs(haplo4_by_eth)

    # 4. Build diplotype freq dfs
    print("Building diplotype frequency tables...")
    diplo5_by_eth = {eth: get_diplotype_frequencies(haplo5_by_eth[eth]) for eth in ETHNICITIES}
    diplo4_by_eth = {eth: get_diplotype_frequencies(haplo4_by_eth[eth]) for eth in ETHNICITIES}
    # The combined pool is kept at HAPLOTYPE level. At a 1e-6 floor it is
    # ~25,000 haplotypes (~300M diplotype pairs) and expanding it would
    # exhaust memory; the cross-ethnic model needs only per-haplotype
    # frequencies (see _align_cross), and the Combined-as-own-population
    # curve uses the numeric vector path below.

    # 5. Run coverage sweeps
    curve_rows = []
    target_rows = []

    for match_level, diplo_by_eth, combined_haplo in [
        ('8of8',  diplo4_by_eth, combined4),
        ('10of10', diplo5_by_eth, combined5),
    ]:
        print(f"\n=== {match_level} ===")

        # Per-ethnicity
        for eth in ETHNICITIES:
            patient_diplo = diplo_by_eth[eth]
            donor_same = diplo_by_eth[eth]
            donor_cross = combined_haplo

            # same_ethnicity
            cov_same = run_sweep(patient_diplo, donor_same, N_SWEEP, cross=False)
            for n, cov in zip(N_SWEEP, cov_same):
                curve_rows.append({
                    'N': n, 'coverage': cov, 'match_level': match_level,
                    'ethnicity': eth, 'model_variant': 'same_ethnicity',
                })
            for thr in THRESHOLDS:
                rs = find_registry_size(patient_diplo, thr)
                target_rows.append({
                    'match_level': match_level, 'ethnicity': eth,
                    'model_variant': 'same_ethnicity',
                    'target_coverage': thr, 'registry_size': rs,
                })

            # cross_ethnic
            cov_cross = run_sweep(patient_diplo, donor_cross, N_SWEEP, cross=True)
            for n, cov in zip(N_SWEEP, cov_cross):
                curve_rows.append({
                    'N': n, 'coverage': cov, 'match_level': match_level,
                    'ethnicity': eth, 'model_variant': 'cross_ethnic',
                })
            for thr in THRESHOLDS:
                rs = find_registry_size_cross(patient_diplo, donor_cross, thr)
                target_rows.append({
                    'match_level': match_level, 'ethnicity': eth,
                    'model_variant': 'cross_ethnic',
                    'target_coverage': thr, 'registry_size': rs,
                })

            print(f"  Processed {eth} {match_level}")

        # Combined (same_ethnicity only — donor=patient=combined).
        # Numeric vector path: needs only diplotype frequencies, so the
        # ~300M-pair combined frame is never materialised with string columns.
        comb_vec = diplotype_freq_vector(combined_haplo['frequency'].to_numpy())
        cov_comb = np.array([
            float(np.clip((comb_vec * (1.0 - np.power(1.0 - comb_vec, float(n)))).sum(), 0.0, 1.0))
            for n in N_SWEEP])
        for n, cov in zip(N_SWEEP, cov_comb):
            curve_rows.append({
                'N': n, 'coverage': cov, 'match_level': match_level,
                'ethnicity': 'Combined', 'model_variant': 'same_ethnicity',
            })
        for thr in THRESHOLDS:
            rs = find_registry_size_vec(comb_vec, thr)
            target_rows.append({
                'match_level': match_level, 'ethnicity': 'Combined',
                'model_variant': 'same_ethnicity',
                'target_coverage': thr, 'registry_size': rs,
            })
        print(f"  Processed Combined {match_level}")

    # 6. Save CSVs
    curves_df = pd.DataFrame(curve_rows)
    curves_df.to_csv(os.path.join(DATA_DIR, 'coverage_curves.csv'), index=False)
    print(f"\nSaved coverage_curves.csv: {len(curves_df)} rows")

    targets_df = pd.DataFrame(target_rows)
    targets_df.to_csv(os.path.join(DATA_DIR, 'registry_size_targets.csv'), index=False)
    print(f"Saved registry_size_targets.csv: {len(targets_df)} rows")

    # 7. Generate figures
    print("\nGenerating figures...")
    for match_level in ['8of8', '10of10']:
        ml_curves = [r for r in curve_rows if r['match_level'] == match_level]
        # Group by ethnicity+model_variant into arrays
        grouped = {}
        for r in ml_curves:
            key = (r['ethnicity'], r['model_variant'])
            if key not in grouped:
                grouped[key] = {'N': [], 'coverage': []}
            grouped[key]['N'].append(r['N'])
            grouped[key]['coverage'].append(r['coverage'])

        fig_rows = []
        for (eth, mv), data in grouped.items():
            fig_rows.append({
                'ethnicity': eth, 'model_variant': mv,
                'N': np.array(data['N']), 'coverage': np.array(data['coverage']),
            })

        fig_path = os.path.join(FIG_DIR, f'coverage_curves_{match_level}.png')
        make_figure(fig_rows, match_level, fig_path)

    # 8. Print verification: Chinese 10of10 same_ethnicity
    print("\n=== Verification: Chinese 10of10 same_ethnicity ===")
    check = targets_df[
        (targets_df['match_level'] == '10of10') &
        (targets_df['ethnicity'] == 'Chinese') &
        (targets_df['model_variant'] == 'same_ethnicity')
    ][['target_coverage', 'registry_size']]
    print(check.to_string(index=False))


if __name__ == '__main__':
    main()
