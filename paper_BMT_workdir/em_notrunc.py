"""Confirmatory: re-run the project's own EM with and without the 0.1% truncation."""
import sys, os, pandas as pd
HERE = '/data/alvin/HLA/analysis'
sys.path.insert(0, HERE)
os.chdir(HERE)
from hwe_test import run_em_haplotypes
from registry_model import get_diplotype_frequencies, find_registry_size, compute_coverage

h = pd.read_csv('data/hla_clean.csv')
h = h[h.ethnicity.isin(['Chinese', 'Others'])]
out = []
for thr in (0.001, 0.0):
    em = run_em_haplotypes(h, freq_threshold=thr, cap=5000)
    for eth in ['Chinese', 'Others']:
        s = em[em.ethnicity == eth][['haplotype', 'frequency']].copy()
        mass = s.frequency.sum()
        s['frequency'] = s.frequency / mass
        dip = get_diplotype_frequencies(s)
        rec = dict(ethnicity=eth, threshold=thr, n_haps=len(s),
                   mass_retained=round(float(mass), 4),
                   N95=find_registry_size(dip, 0.95),
                   cov_at_50k=round(100 * compute_coverage(dip, 50000), 1))
        out.append(rec)
        print(rec, flush=True)
pd.DataFrame(out).to_csv('/data/alvin/HLA/paper_BMT_workdir/em_notrunc.csv', index=False)
print("DONE")
