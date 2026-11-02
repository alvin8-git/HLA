"""Coverage within the validated extrapolation range (2x, 3x current size), and
model N* from the live 1e-6 run, for the review's recommendation table.
Usage: python publishable_range.py <haplo_freqs_em.csv> <dedup_coverage_targets.csv> <live registry_size_targets.csv> <snapshot registry_size_targets.csv>
"""
import sys
import numpy as np, pandas as pd
h = pd.read_csv(sys.argv[1]); t = pd.read_csv(sys.argv[2]); live = pd.read_csv(sys.argv[3]); snap = pd.read_csv(sys.argv[4])


def coverage(f, N):
    f = np.sort(f)[::-1]; tot = 0.0
    for i in range(len(f)):
        d = 2.0 * f[i] * f[i:]; d[0] = f[i] ** 2
        tot += np.sum(d * -np.expm1(N * np.log1p(-d)))
    return tot


def cj_cov(n, f1, f2, N):
    A = (n - 1) * f1 / ((n - 1) * f1 + 2 * f2); return 1 - f1 / n * A ** (N - n + 1)


r = t[t.scenario == 'near'].set_index('ethnicity')
print("coverage at k x current (near-dedup) size: genotype Chao-Jost | 1e-6 haplotype model")
for eth in ['Chinese', 'Malay', 'Indian', 'Others']:
    n, f1, f2 = int(r.loc[eth, 'n_bmdp']), int(r.loc[eth, 'f1']), int(r.loc[eth, 'f2'])
    f = h[h.ethnicity == eth].frequency.to_numpy()
    s = "  ".join(f"{k}x: {cj_cov(n, f1, f2, k*n):.3f} | {coverage(f, k*n):.3f}" for k in (1, 2, 3))
    print(f"{eth:8s} n={n:6d}  {s}")
sel = lambda d: d[(d.model_variant == 'same_ethnicity')].pivot_table(index=['match_level', 'ethnicity'], columns='target_coverage', values='registry_size')
print("\nlive 1e-6 model N* (same-ethnicity):"); print(sel(live).round(0).astype('Int64').to_string())
print("\nsnapshot 1e-3 model N*:"); print(sel(snap).round(0).astype('Int64').to_string())
