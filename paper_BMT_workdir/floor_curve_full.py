"""Floor sweep at the FULL sample (cap=50000), superseding the cap=5000 sweep.

The earlier sweep ran the EM on 5,000 individuals per group. 15_em_convergence
later showed that cap inflates N* by 264% at a 1e-4 floor, because a 5,000-
individual EM retains spurious phase-ambiguity haplotypes. Any conclusion about
whether the rare tail is real signal therefore has to be re-derived here.

Key scale reference: a singleton haplotype in a sample of n individuals (2n
chromosomes) has frequency 1/(2n). For Chinese n=44,400 that is 1.1e-5, so the
floors below bracket the singleton threshold rather than sitting above it.
"""
import sys, os, gc, numpy as np, pandas as pd

HERE = '/data/alvin/HLA/analysis'
sys.path.insert(0, HERE)
os.chdir(HERE)
from hwe_test import run_em_haplotypes
from registry_model import diplotype_freq_vector, find_registry_size_vec

FLOORS = [0.0, 1e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3]
ETHS = ['Chinese', 'Others']

h = pd.read_csv('data/hla_clean.csv')
h = h[h.ethnicity.isin(ETHS)]

print("running unfloored EM at full sample (cap=50000)...", flush=True)
em = run_em_haplotypes(h, freq_threshold=0.0, cap=50000)
print("  done:", em.groupby('ethnicity').size().to_dict(), flush=True)

n_indiv = h.groupby('ethnicity')['sample_id'].nunique().to_dict()
rows = []
for eth in ETHS:
    full = em[em.ethnicity == eth]['frequency'].to_numpy(dtype=np.float64)
    full = full / full.sum()
    singleton = 1.0 / (2 * n_indiv[eth])
    print(f"\n{eth}: n={n_indiv[eth]:,}  singleton freq={singleton:.2e}  "
          f"unfloored haplotypes={len(full):,}", flush=True)
    for fl in FLOORS:
        kept = full[full >= fl] if fl > 0 else full
        mass = float(kept.sum())
        vec = diplotype_freq_vector(kept)
        cov50 = float(np.clip((vec * (1 - np.power(1 - vec, 5e4))).sum(), 0, 1))
        cov500 = float(np.clip((vec * (1 - np.power(1 - vec, 5e5))).sum(), 0, 1))
        n90 = find_registry_size_vec(vec, 0.90)
        n95 = find_registry_size_vec(vec, 0.95)
        rec = dict(ethnicity=eth, n_individuals=n_indiv[eth], floor=fl,
                   n_haps=int(len(kept)), mass_retained=round(mass, 4),
                   n_diplo=int(len(vec)),
                   cov_50k=round(100 * cov50, 1), cov_500k=round(100 * cov500, 1),
                   N90=n90, N95=n95)
        rows.append(rec)
        print(f"  floor={fl:<8g} haps={len(kept):>7,} mass={mass:6.1%} "
              f"cov@50k={100*cov50:5.1f}% N90={n90:>14,} N95={n95:>14,}", flush=True)
        del vec; gc.collect()

df = pd.DataFrame(rows)
df.to_csv('/data/alvin/HLA/paper_BMT_workdir/floor_curve_full.csv', index=False)
print("\nDONE -> floor_curve_full.csv")
