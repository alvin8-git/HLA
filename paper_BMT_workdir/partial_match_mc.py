"""Monte-Carlo partial-match coverage at the 1e-4 floor.

The exact algorithm in 06_partial_match_plots.py is O(N_diplo^2). At a 1e-3
floor N_diplo ~ 1e4 and that is fine; at 1e-4 it is ~2.8e6 and the exact form
needs ~8e12 diplotype comparisons.

Here patients are sampled in proportion to diplotype frequency and compared
against the FULL donor diplotype set, which is unbiased for
    C(N) = E_patient[ 1 - (1 - p_match)^N ]
with Monte-Carlo error O(1/sqrt(S)) rather than any truncation bias.
"""
import sys, os, importlib.util, numpy as np, pandas as pd

HERE = '/data/alvin/HLA/analysis'
sys.path.insert(0, HERE)
os.chdir(HERE)
spec = importlib.util.spec_from_file_location("pm", "06_partial_match_plots.py")
pm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pm)

S = 3000            # sampled patient diplotypes
CHUNK = 4000        # donor diplotypes per chunk
SEED = 42
THRESHOLDS = (0.75, 0.90, 0.95)
ETHS = sys.argv[1:] or ['Chinese']

haplo = pd.read_csv('data/haplo_freqs_em.csv')
rng = np.random.default_rng(SEED)
out = []

for eth in ETHS:
    haps = pm.parse_haplotypes(haplo[haplo.ethnicity == eth], pm.LOCI_5)
    alleles = [h for h, f in haps]
    freqs = np.array([f for h, f in haps], dtype=np.float64)
    freqs = freqs / freqs.sum()
    vocab = sorted({a for h in alleles for a in h})
    enc = {a: np.int16(i) for i, a in enumerate(vocab)}
    arr = np.array([[enc[a] for a in h] for h in alleles], dtype=np.int16)

    ii, jj = np.triu_indices(len(alleles))
    f_dip = 2.0 * freqs[ii] * freqs[jj]
    same = ii == jj
    f_dip[same] = freqs[ii[same]] ** 2
    f_dip = f_dip / f_dip.sum()
    print(f"{eth}: {len(alleles)} haplotypes -> {len(f_dip):,} diplotypes", flush=True)

    pat_idx = rng.choice(len(f_dip), size=S, replace=True, p=f_dip)
    pa, pb = arr[ii[pat_idx]], arr[jj[pat_idx]]

    acc = {m: np.zeros(S) for m in (10, 9, 8)}
    for start in range(0, len(f_dip), CHUNK):
        end = min(start + CHUNK, len(f_dip))
        cm = pm.allele_match_count_matrix(pa, pb, arr[ii[start:end]], arr[jj[start:end]])
        w = f_dip[start:end]
        for m in acc:
            acc[m] += (cm >= m) @ w
        if (start // CHUNK) % 200 == 0:
            print(f"  {eth} {100*end/len(f_dip):5.1f}%", flush=True)

    for m in (10, 9, 8):
        p = np.clip(acc[m], 0.0, 1.0)
        for thr in THRESHOLDS:
            lo, hi = 1.0, 1e10
            for _ in range(60):
                mid = (lo * hi) ** 0.5
                cov = float(np.mean(1.0 - np.power(1.0 - p, mid)))
                lo, hi = (mid, hi) if cov < thr else (lo, mid)
            out.append(dict(ethnicity=eth, min_match=f"{m}/10", coverage=thr,
                            N_star=int(round(hi))))
        cov50k = float(np.mean(1.0 - np.power(1.0 - p, 50_000)))
        out.append(dict(ethnicity=eth, min_match=f"{m}/10", coverage='cov@50k',
                        N_star=round(100 * cov50k, 1)))
        print(f"  {eth} >={m}/10: coverage at 50k = {100*cov50k:.1f}%", flush=True)

pd.DataFrame(out).to_csv('/data/alvin/HLA/paper_BMT_workdir/partial_match_mc.csv', index=False)
print("DONE")
