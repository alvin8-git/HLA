"""Regenerate Figures 3 and 4 (partial-match coverage) at the corrected floor.

06_partial_match_plots.py is O(N_diplo^2) and cannot run once a 1e-4 floor
pushes the diplotype count into the millions. Patients are therefore sampled in
proportion to diplotype frequency and evaluated against the FULL donor diplotype
set, which is unbiased for coverage with Monte-Carlo error ~1/sqrt(S) and no
truncation bias.

Outputs (overwriting the stale 1e-3-era figures):
  analysis/figures/partial_match_10locus.png
  analysis/figures/partial_match_8locus.png
  paper_BMT_workdir/partial_match_mc_all.csv
"""
import sys, os, importlib.util, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = '/data/alvin/HLA/analysis'
sys.path.insert(0, HERE)
os.chdir(HERE)
spec = importlib.util.spec_from_file_location("pm", "06_partial_match_plots.py")
pm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pm)

S, CHUNK, SEED = 3000, 4000, 42
ETHS = ['Chinese', 'Malay', 'Indian', 'Others']
N_GRID = np.logspace(3, 9, 60)
FRAMEWORKS = [
    ('10-locus (HLA-A, B, C, DRB1, DQB1) matching', pm.LOCI_5,
     [(10, '10/10', 'green'), (9, '≥9/10', 'blue'), (8, '≥8/10', 'red')],
     'partial_match_10locus.png'),
    ('8-locus (HLA-A, B, C, DRB1) matching', pm.LOCI_4,
     [(8, '8/8', 'green'), (7, '≥7/8', 'blue'), (6, '≥6/8', 'red')],
     'partial_match_8locus.png'),
]

haplo = pd.read_csv('data/haplo_freqs_em.csv')
rng = np.random.default_rng(SEED)
rows = []

for fw_name, loci, levels, out_png in FRAMEWORKS:
    print(f"\n=== {fw_name} ===", flush=True)
    curves = {}
    for eth in ETHS:
        haps = pm.parse_haplotypes(haplo[haplo.ethnicity == eth], loci)
        alleles = [h for h, f in haps]
        freqs = np.array([f for h, f in haps], dtype=np.float64)
        freqs /= freqs.sum()
        vocab = sorted({a for h in alleles for a in h})
        enc = {a: np.int16(i) for i, a in enumerate(vocab)}
        arr = np.array([[enc[a] for a in h] for h in alleles], dtype=np.int16)

        ii, jj = np.triu_indices(len(alleles))
        f_dip = 2.0 * freqs[ii] * freqs[jj]
        same = ii == jj
        f_dip[same] = freqs[ii[same]] ** 2
        f_dip /= f_dip.sum()

        pat = rng.choice(len(f_dip), size=S, replace=True, p=f_dip)
        pa, pb = arr[ii[pat]], arr[jj[pat]]
        acc = {m: np.zeros(S) for m, _, _ in levels}
        for st in range(0, len(f_dip), CHUNK):
            en = min(st + CHUNK, len(f_dip))
            cm = pm.allele_match_count_matrix(pa, pb, arr[ii[st:en]], arr[jj[st:en]])
            w = f_dip[st:en]
            for m in acc:
                acc[m] += (cm >= m) @ w
        curves[eth] = {}
        for m, label, _ in levels:
            p = np.clip(acc[m], 0.0, 1.0)
            cov = np.array([float(np.mean(1.0 - np.power(1.0 - p, n))) for n in N_GRID])
            curves[eth][m] = cov
            for thr in (0.75, 0.90, 0.95):
                lo, hi = 1.0, 1e10
                for _ in range(60):
                    mid = (lo * hi) ** 0.5
                    c = float(np.mean(1.0 - np.power(1.0 - p, mid)))
                    lo, hi = (mid, hi) if c < thr else (lo, mid)
                rows.append(dict(framework=fw_name.split()[0], ethnicity=eth,
                                 min_match=label, coverage=thr, N_star=int(round(hi))))
            rows.append(dict(framework=fw_name.split()[0], ethnicity=eth,
                             min_match=label, coverage='cov@50k',
                             N_star=round(100 * float(np.mean(1 - np.power(1 - p, 5e4))), 1)))
        print(f"  {eth}: {len(alleles)} haps, {len(f_dip):,} diplotypes", flush=True)

    fig, axes = plt.subplots(1, 4, figsize=(19, 4.2), sharey=True, facecolor='white')
    for ax, eth in zip(axes, ETHS):
        for m, label, colour in levels:
            ax.plot(N_GRID, 100 * curves[eth][m], color=colour, lw=1.9, label=label)
        for y in (75, 90, 95):
            ax.axhline(y, color='grey', ls='--', lw=0.7, alpha=0.6)
        ax.axvline(5e4, color='black', ls=':', lw=1.1)
        ax.set_xscale('log'); ax.set_xlim(1e3, 1e9); ax.set_ylim(0, 100)
        ax.set_title(eth, fontweight='bold'); ax.set_xlabel('Registry size (N donors)')
        ax.set_facecolor('white')
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel('Patients with ≥1 matched donor (%)')
    axes[0].legend(fontsize=9, loc='upper left', frameon=False)
    axes[0].annotate('N = 50,000', xy=(5e4, 4), fontsize=8, rotation=90, va='bottom')
    fig.suptitle(f'{fw_name} — Monte-Carlo estimate, 1e-4 floor, full-sample EM',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    out = os.path.join('figures', out_png)
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: {out}", flush=True)

pd.DataFrame(rows).to_csv('/data/alvin/HLA/paper_BMT_workdir/partial_match_mc_all.csv', index=False)
print("\nDONE")
