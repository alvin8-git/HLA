"""1e-6 haplotype-HWE model (unconditional, no pooling) at the deduplicated BMDP sizes,
against the SCBB probe truth. Also times one full-tail coverage evaluation.
Usage: python model_vs_probe.py <haplo_freqs_em.csv> <dedup_coverage_targets.csv>
"""
import sys, time
import numpy as np, pandas as pd
h = pd.read_csv(sys.argv[1]); t = pd.read_csv(sys.argv[2])


def coverage(f, N):
    f = np.sort(f)[::-1]; tot = 0.0
    for i in range(len(f)):
        d = 2.0 * f[i] * f[i:]; d[0] = f[i] ** 2
        tot += np.sum(d * -np.expm1(N * np.log1p(-d)))
    return tot


for eth, g in h.groupby('ethnicity'):
    f = g.frequency.to_numpy(); r = t[t.ethnicity == eth].set_index('scenario')
    t0 = time.time(); c_near = coverage(f, r.loc['near', 'n_bmdp']); dt = time.time() - t0
    c_all = coverage(f, r.loc['all', 'n_bmdp'])
    print(f"{eth:8s} K={len(f):5d}  model@near n={int(r.loc['near','n_bmdp'])}: {c_near:.3f}  model@all n={int(r.loc['all','n_bmdp'])}: {c_all:.3f}  "
          f"probe={r.loc['near','probe_cov']:.3f} [{r.loc['near','probe_lo']:.3f},{r.loc['near','probe_hi']:.3f}]  eval {dt:.1f}s")
