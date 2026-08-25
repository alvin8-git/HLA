"""
13_cross_ethnic_sensitivity.py
Cross-ethnic sensitivity analysis: how do registry size targets change under
different ethnic composition assumptions for donors and patients?

Scenarios compared:
  1. SG population weights (Census 2020): Chinese 74.3%, Malay 13.5%, Indian 9.0%, Others 3.2%
  2. BMDP+SCBB registry composition (actual donors): derived from hla_clean.csv source counts
  3. Patient.txt composition: ethnic breakdown of actual patients seeking donors
  4. Worst-case minority-focus: Indian 50%, Malay 30%, Others 20%, Chinese 0%

For each scenario, donor and patient pools can be set independently.
Main analysis: same-composition (donors = patients = scenario weights).
Cross composition: patient weights fixed at SG pop; vary donor weights.

Outputs:
  data/cross_ethnic_sensitivity.csv    — N* per scenario × threshold
  figures/cross_ethnic_sensitivity.png — bar chart comparison
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from registry_model import get_diplotype_frequencies, find_registry_size

DATA_DIR = os.path.join(HERE, 'data')
FIG_DIR  = os.path.join(HERE, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

ETHNICITIES = ['Chinese', 'Malay', 'Indian', 'Others']
THRESHOLDS  = [0.75, 0.85, 0.90, 0.95]

ETH_MAP = {'BMDP_OUT': 'Chinese',   # predominantly Chinese
            'SCBB_OUT': 'Chinese',   # check hla_clean for actual breakdown
            'HSA-Donor': 'Others',
            'HSA-Patient': 'Others'}


# ---------------------------------------------------------------------------
# Derive registry composition from hla_clean.csv
# ---------------------------------------------------------------------------

def get_registry_composition(hla_df: pd.DataFrame) -> dict:
    """Compute ethnic proportions of BMDP+SCBB donors (excluding HSA)."""
    donors = hla_df[hla_df['source'].isin(['BMDP_OUT', 'SCBB_OUT'])]
    counts = donors.groupby('ethnicity')['sample_id'].nunique()
    total = counts.sum()
    weights = {eth: round(counts.get(eth, 0) / total, 4) for eth in ETHNICITIES}
    return weights


def get_patient_composition(patient_path: str) -> dict:
    """Compute ethnic proportions from Patient.txt."""
    df = pd.read_csv(patient_path, sep='\t', header=None,
                     names=['eth_code', 'A', 'B', 'C', 'DRB1', 'DQB1'])
    eth_map = {'C': 'Chinese', 'M': 'Malay', 'I': 'Indian', 'O': 'Others'}
    df['ethnicity'] = df['eth_code'].map(eth_map)
    counts = df['ethnicity'].value_counts()
    total = counts.sum()
    weights = {eth: round(counts.get(eth, 0) / total, 4) for eth in ETHNICITIES}
    return weights


# ---------------------------------------------------------------------------
# Coverage computation
# ---------------------------------------------------------------------------

def compute_same_ethnic_weighted(weights: dict, haplo_dfs: dict,
                                  threshold: float) -> int:
    """
    N* assuming same-ethnicity matching only.
    N* = weighted average: Σ_eth w_eth × N*_eth
    (conservative approximation for cross-ethnic).
    """
    total_n = 0.0
    for eth in ETHNICITIES:
        w = weights.get(eth, 0.0)
        if w < 1e-6:
            continue
        df = haplo_dfs.get(eth)
        if df is None or df.empty:
            continue
        diplo = get_diplotype_frequencies(df)
        n_eth = find_registry_size(diplo, threshold)
        total_n += w * n_eth
    return int(round(total_n))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Loading EM haplotype frequencies...")
    em_df = pd.read_csv(os.path.join(DATA_DIR, 'haplo_freqs_em.csv'))

    haplo_dfs = {}
    for eth in ETHNICITIES:
        sub = em_df[em_df['ethnicity'] == eth][['haplotype', 'frequency']].copy()
        sub['frequency'] = sub['frequency'] / sub['frequency'].sum()
        haplo_dfs[eth] = sub

    print("Loading hla_clean.csv for registry composition...")
    hla_df = pd.read_csv(os.path.join(DATA_DIR, 'hla_clean.csv'))

    registry_weights = get_registry_composition(hla_df)
    print(f"  Registry (BMDP+SCBB) weights: {registry_weights}")

    root = os.path.dirname(HERE)
    patient_weights = get_patient_composition(os.path.join(root, 'Patient.txt'))
    print(f"  Patient.txt weights:          {patient_weights}")

    # Define scenarios
    scenarios = {
        'SG population\n(current model)': {
            'Chinese': 0.743, 'Malay': 0.135, 'Indian': 0.090, 'Others': 0.032},
        'BMDP+SCBB\nregistry composition': registry_weights,
        'Patient.txt\ncomposition': patient_weights,
        'Minority-focus\n(Indian+Malay+Others)': {
            'Chinese': 0.00, 'Malay': 0.40, 'Indian': 0.40, 'Others': 0.20},
    }

    print("\nComputing N* per scenario (same-ethnicity weighted average)...")
    rows = []
    for scenario_name, weights in scenarios.items():
        print(f"\n  Scenario: {scenario_name.replace(chr(10), ' ')}")
        print(f"    Weights: {weights}")
        for t in THRESHOLDS:
            n = compute_same_ethnic_weighted(weights, haplo_dfs, t)
            print(f"    N* at {int(t*100)}%: {n:,}")
            rows.append({
                'scenario': scenario_name.replace('\n', ' '),
                'target_coverage': t,
                'registry_size': n,
                **{f'weight_{eth}': weights.get(eth, 0) for eth in ETHNICITIES},
            })

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(DATA_DIR, 'cross_ethnic_sensitivity.csv'), index=False)
    print(f"\nSaved: cross_ethnic_sensitivity.csv")

    # -----------------------------------------------------------------------
    # Figure: grouped bar chart, 4 thresholds × 4 scenarios
    # -----------------------------------------------------------------------
    threshold_labels = {0.75: '75%', 0.85: '85%', 0.90: '90%', 0.95: '95%'}
    scenario_names   = [s.replace('\n', ' ') for s in scenarios.keys()]
    scenario_colors  = ['#2c7bb6', '#d7191c', '#1a9641', '#fdae61']

    x = np.arange(len(THRESHOLDS))
    bar_width = 0.18
    n_scen = len(scenario_names)

    fig, ax = plt.subplots(figsize=(11, 5), facecolor='white')
    ax.set_facecolor('white')

    for k, sname in enumerate(scenario_names):
        sub = df[df['scenario'] == sname]
        vals = [sub[sub['target_coverage'] == t]['registry_size'].values[0]
                for t in THRESHOLDS]
        offset = (k - (n_scen - 1) / 2) * bar_width
        bars = ax.bar(x + offset, vals, bar_width,
                      label=sname, color=scenario_colors[k], alpha=0.85,
                      edgecolor='white')
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 200,
                    f'{int(v):,}', ha='center', va='bottom',
                    fontsize=6.5, rotation=90)

    ax.set_xticks(x)
    ax.set_xticklabels([threshold_labels[t] for t in THRESHOLDS], fontsize=11)
    ax.set_xlabel('Coverage target', fontsize=11)
    ax.set_ylabel('Estimated registry size (donors)', fontsize=11)
    ax.set_title('Cross-ethnic sensitivity analysis: registry size targets\nby ethnic composition scenario (same-ethnicity weighted average)',
                 fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_ylim(0, ax.get_ylim()[1] * 1.25)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{int(v):,}'))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for y in [10000, 20000, 30000, 40000, 50000]:
        ax.axhline(y, color='lightgrey', linestyle='--', linewidth=0.7, zorder=0)

    plt.tight_layout()
    out_path = os.path.join(FIG_DIR, 'cross_ethnic_sensitivity.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {out_path}')

    print('\n=== N* at 95% coverage by scenario ===')
    mask = df['target_coverage'] == 0.95
    print(df[mask][['scenario', 'registry_size']].to_string(index=False))
    print('\nDone.')


if __name__ == '__main__':
    main()
