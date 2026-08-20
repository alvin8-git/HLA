#!/usr/bin/env bash
# Overnight chain v4 (adds 14 + 15, the two figure scripts omitted from v3): uncapped EM at a 1e-6 floor -> downstream -> B=1000 bootstrap -> v2.16.
#
# Supersedes the earlier chain, which was started before 15_em_convergence.py
# revealed that the 5,000-individual EM cap inflates the Chinese N* by 264%
# at a 1e-4 floor (11,487,962 capped vs 3,153,571 at the full 45,018 sample).
# The cap is now 50,000 and no longer binds for any CMIO group, so the whole
# chain must be recomputed from the EM step down.
#
# Launched detached; survives session end. Progress: paper_BMT_workdir/overnight.log
set -u
cd /data/alvin/HLA
# Pin the interpreter. A login shell (screen, cron, at) puts /usr/bin/python
# (2.7) on PATH ahead of conda, which fails on the first non-ASCII docstring.
PY="${PY:-/home/alvin/miniconda3/bin/python}"
[ -x "$PY" ] || { echo "no interpreter at $PY" >&2; exit 127; }
WD=paper_BMT_workdir
LOG=$WD/overnight.log
say() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >> "$LOG"; }

step() {           # step <label> <logfile> <dir> <cmd...>
  local label=$1 logf=$2 dir=$3; shift 3
  say "-> $label"
  ( cd "$dir" && "$@" ) > "$WD/$logf" 2>&1
  local rc=$?
  if [ $rc -ne 0 ]; then
    say "ERROR: $label exited $rc"
    say "       $(tail -3 $WD/$logf | tr '\n' ' ')"
    say "=== chain aborted ==="
    exit $rc
  fi
  say "   $label ok"
}

rm -f $WD/OVERNIGHT_DONE
say "=== overnight chain v3 started (uncapped EM, floor 1e-6, B=1000) ==="

# 1. EM haplotype frequencies, cap=50000, floor=1e-6 (below every group's 1/(2n))
step "03 EM haplotypes (uncapped, 1e-6)" run03_1e6.log analysis "$PY" -u 03_hwe_test.py
say "   $("$PY" - <<'PY'
import pandas as pd
h=pd.read_csv('analysis/data/haplo_freqs_em.csv')
g=h.groupby('ethnicity').agg(n=('frequency','size'),mass=('frequency','sum'))
print('; '.join(f"{e}: {r.n} haps, {r.mass:.1%} mass" for e,r in g.iterrows()))
PY
)"

# 2. Core model + downstream analyses that read the haplotype table
step "04 registry model"        run04_uncapped.log analysis "$PY" -u 04_registry_model.py
step "13 cross-ethnic sensitivity" run13_uncapped.log analysis "$PY" -u 13_cross_ethnic_sensitivity.py
step "16 smoothing sensitivity" run16_uncapped.log analysis "$PY" -u 16_smoothing_sensitivity.py
step "10 LD report"              run10_1e6.log analysis "$PY" -u 10_ld_report.py
step "11 Others stratification" run11_1e6.log analysis "$PY" -u 11_others_stratification.py
step "07 GENE[RATE] validation" run07_uncapped.log analysis "$PY" -u 07_validate_em.py
step "12 match validation"      run12_uncapped.log analysis "$PY" -u 12_match_validation.py
step "08 section plots"         run08_uncapped.log analysis "$PY" -u 08_section_plots.py
step "14 pipeline flowchart"   run14.log       analysis "$PY" -u 14_pipeline_flowchart.py
step "15 EM convergence"       run15.log       analysis "$PY" -u 15_em_convergence.py
step "partial-match figures+CSV" pm_figs_1e6.log . "$PY" paper_BMT_workdir/partial_match_figs.py

# 3. Full bootstrap
say "-> 09 bootstrap B=1000 (expect ~3-4h)"
( cd analysis && HLA_BOOTSTRAP_B=1000 "$PY" -u 09_bootstrap_ci.py ) > $WD/run09_B1000.log 2>&1
RC=$?
if [ $RC -ne 0 ]; then
  say "ERROR: bootstrap exited $RC"; say "       $(tail -3 $WD/run09_B1000.log | tr '\n' ' ')"
  say "=== chain aborted ==="; exit $RC
fi
say "   bootstrap ok"

# 4. Manuscript
step "build v2.16" run_build_v216.log . "$PY" build_report.py

# 5. Sanity report
"$PY" - >> "$LOG" 2>&1 <<'PY'
import pandas as pd
ci = pd.read_csv('analysis/data/registry_size_ci.csv')
s = ci[(ci.match_level=='10of10') & (ci.target_coverage==0.95)]
print("    FINAL 10/10 @95% (bootstrap median [95% CI]):")
for r in s.itertuples():
    print(f"      {r.ethnicity:8s} {r.registry_size:>14,}  [{r.ci_lo:,} - {r.ci_hi:,}]")
t = pd.read_csv('analysis/data/registry_size_targets.csv')
x = t[(t.match_level=='10of10') & (t.target_coverage==0.95) &
      (t.model_variant=='same_ethnicity')]
print("    EM point estimates (same-ethnicity):")
for r in x.itertuples():
    print(f"      {r.ethnicity:8s} {int(r.registry_size):>14,}")
PY
say "=== overnight chain finished OK ==="
touch $WD/OVERNIGHT_DONE
