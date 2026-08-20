#!/usr/bin/env bash
# Run this once the overnight chain reports done. It gates on artifact
# freshness, prints the headline numbers for you to eyeball, then commits.
#
#   cd /data/alvin/HLA && ./finish_and_commit.sh
#
# Aborts rather than committing if anything is stale. Nothing here is
# irreversible until the final `git commit`, which is the last line.
set -u
cd /data/alvin/HLA
PY="${PY:-/home/alvin/miniconda3/bin/python}"

echo "== 1. chain status =="
if [ -f paper_BMT_workdir/OVERNIGHT_DONE ]; then
  echo "   OVERNIGHT_DONE present"
else
  echo "   NOT FINISHED — last log lines:"; tail -4 paper_BMT_workdir/overnight.log
  echo "   (attach with: screen -r hla)"; exit 1
fi

echo
echo "== 2. freshness gate =="
"$PY" check_freshness.py --quiet || {
  echo "   STALE ARTIFACTS — not committing. Re-run the offending scripts above."; exit 1; }

echo
echo "== 3. headline numbers (sanity-check these before the commit lands) =="
"$PY" - <<'EOF'
import pandas as pd
t = pd.read_csv('analysis/data/registry_size_targets.csv')
s = t[(t.match_level=='10of10') & (t.model_variant=='same_ethnicity') &
      (t.target_coverage==0.95)]
print("  10/10 @95% same-ethnicity:")
for r in s.itertuples():
    print(f"    {r.ethnicity:9s} {int(r.registry_size):>14,}")
try:
    ci = pd.read_csv('analysis/data/registry_size_ci.csv')
    c = ci[(ci.match_level=='10of10') & (ci.target_coverage==0.95)]
    print("  bootstrap median [95% CI]:")
    for r in c.itertuples():
        print(f"    {r.ethnicity:9s} {r.registry_size:>14,}  [{r.ci_lo:,} - {r.ci_hi:,}]")
except Exception as e:
    print("  (no CI file:", e, ")")
h = pd.read_csv('analysis/data/haplo_freqs_em.csv')
g = h.groupby('ethnicity').agg(n=('frequency','size'), mass=('frequency','sum'))
print("  haplotypes / mass retained:",
      "; ".join(f"{e} {int(r.n):,} ({r.mass:.1%})" for e, r in g.iterrows()))
EOF

echo
echo "== 4. staging =="
git add -A analysis/ build_report.py build_paper_docx.py check_freshness.py \
          run_overnight.sh finish_and_commit.sh VERSION.md .gitignore CLAUDE.md \
          AUDIT_v2.15_findings.md HLA_Registry_Size_CMIO_BoneMarrowTransplantation.md \
          paper_BMT_workdir/
# excluded deliberately: cartoons + slide deck (Alvin: "do not commit the cartoons"),
# *.docx (gitignored, regenerated per release), raw .xlsx/.pdf source data
# (untracked by repo convention), graphify-out/, .claude/, .omc/, BACKUP_*.csv
git reset -q paper_BMT_workdir/*.log 2>/dev/null
git status --short | head -40
echo "  (files staged: $(git diff --cached --name-only | wc -l))"

echo
echo "== 5. commit =="
git commit -q -F - <<'MSG'
fix: rare-haplotype floor must sit below 1/(2n) — retracts the v2.15 registry targets

The ~40,000-45,000 donors-per-CMIO-group figure reported through v2.15 is an
artefact of preprocessing, not a property of the population. Two parameters
that were never reported dominate the result.

A frequency floor is harmless below the singleton frequency 1/(2n) and
destructive above it. For Chinese donors (n=45,754, singleton 1.09e-5),
discarding the 225k haplotypes below 1e-6 — 96% of all distinct haplotypes —
costs 0.03% of frequency mass and leaves N* unchanged; raising the floor to
1e-3 collapses N* 2,098-fold. The Others group (singleton 1.27e-4) reproduces
the rule with its threshold displaced accordingly. The sub-singleton tail is
EM phase-ambiguity noise and is provably inert, so the damage is done to
haplotypes the sample genuinely resolves.

The EM input cap interacts with the floor: a 5,000-individual cap costing ~8%
at 1e-3 inflates by 264% at 1e-4, and its direction of bias is not stable
across floors. Floor and cap cannot be chosen independently.

Because 1/(2n) differs between groups, a single floor biases unequally-sampled
groups unequally — cross-group registry-size rankings are therefore unsafe and
are withdrawn. Within-group comparisons (match relaxation, cross-ethnic
penalty) survive.

Pipeline now runs at freq_threshold=1e-6, cap=50000 (binds for no group),
search ceiling 1e10, retaining 100% of frequency mass in all four groups.

Also fixed:
- cross-ethnic diplotype merge was order-sensitive; pairs are now labelled in
  canonical order (30.4% of patient frequency mass had been scored unmatched)
- bootstrap N_EFF summed to 61,149 against a 59,186 study total
- Table 1's "Weighted Average" row was the pooled-registry model, not an average
- Others fails HWE at all five loci (Wahlund); pooled figures inadmissible
- 10_ld_report reads the haplotype table, so D' is floor-dependent (0.90-0.96)
- 11 and 15 held private FREQ_THRESHOLD copies that ignored the global change
- manuscript prose now interpolates figures from the CSVs instead of hardcoding

Adds check_freshness.py, which fails the build if any derived artifact predates
the haplotype table it comes from, and run_overnight.sh to drive the pipeline
end to end under screen with a pinned interpreter.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
MSG
echo "   committed: $(git log --oneline -1)"
echo
echo "Not committed by design: cartoons/slide deck, *.docx, raw xlsx/pdf, graphify-out/."
echo "Manuscripts are regenerated, not stored: python build_report.py"
