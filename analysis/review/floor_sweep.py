"""Score haplotype-model floors against the observed in-registry 10/10 twin rate.

For each floor: keep EM haplotypes with f >= floor, then predict coverage at the
registry's actual size N = n-1 (leave-one-out), three ways:
  uncond      : no renormalisation, no pooling. A patient carrying a dropped
                haplotype has zero match probability. This is the honest estimand.
  cond        : renormalised to 1, no pooling (coverage among common-haplotype carriers).
  cond_pool99 : renormalised + bottom 1% pooled into one 'other' haplotype, as
                registry_model.get_diplotype_frequencies does in production.
Usage: python floor_sweep.py <haplo_freqs_em.csv> <genotype_coverage_targets.csv> <out.csv>
"""
import sys
import numpy as np, pandas as pd

H, T, OUT = sys.argv[1:4]
h = pd.read_csv(H)
t = pd.read_csv(T)
t = t[t.match_level == '10of10'].set_index('ethnicity')


def coverage(f, N):
    """sum over unordered pairs i<=j of d_ij * [1-(1-d_ij)^N], streamed row by row."""
    f = np.sort(f)[::-1]
    tot = 0.0
    for i in range(len(f)):
        d = 2.0 * f[i] * f[i:]
        d[0] = f[i] ** 2
        tot += np.sum(d * -np.expm1(N * np.log1p(-d)))
    return tot


def pool99(f):
    f = np.sort(f)[::-1]
    cut = min(int((np.cumsum(f) < 0.99).sum()) + 1, len(f))
    res = f[cut:].sum()
    return np.r_[f[:cut], res] if res > 0 else f[:cut]


rows = []
for eth, g in h.groupby('ethnicity'):
    n = int(t.loc[eth, 'n']); obs = t.loc[eth, 'observed_loo_cov']; N = n - 1
    f_all = g.frequency.to_numpy()
    floors = {'1e-6': 1e-6, '1/(4n)': 1 / (8 * n), '1/(2n)': 1 / (2 * n), '2/(2n)': 2 / (2 * n),
              '1e-4': 1e-4, '3e-4': 3e-4, '1e-3': 1e-3}
    for name, fl in floors.items():
        f = f_all[f_all >= fl]
        mass = f.sum()
        fr = f / mass
        r = dict(ethnicity=eth, floor=name, floor_value=fl, K=len(f), mass_retained=mass,
                 observed=obs, uncond=coverage(f, N), cond=coverage(fr, N), cond_pool99=coverage(pool99(fr), N))
        for k in ('uncond', 'cond', 'cond_pool99'):
            r[k + '_err_pts'] = 100 * (r[k] - obs)
        rows.append(r)
        print(f"{eth:8s} {name:7s} K={len(f):5d} mass={mass:.4f} obs={obs:.3f} "
              f"uncond={r['uncond']:.3f} cond={r['cond']:.3f} pool99={r['cond_pool99']:.3f}", flush=True)
df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
print(df.groupby('floor', sort=False)[['uncond_err_pts', 'cond_err_pts', 'cond_pool99_err_pts']]
      .agg(lambda x: np.mean(np.abs(x))).round(2).rename(columns=lambda c: 'MAE_' + c).to_string())
