# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A reproducible Python data-analysis pipeline (no app, no server). Two deliverables:
1. **Verification** — independently re-derives HLA allele/haplotype frequencies from Singapore's BMDP+SCBB bone marrow registry (n=59,186) and audits the 2022 *Blood Cell Therapy* paper (`2022_HLA_BloodCellTherapy.pdf`).
2. **Registry-size model** — minimum donor count for 75–95% patient match probability at 8/8 and 10/10 HLA match, by CMIO ethnicity (Chinese/Malay/Indian/Others) and same-ethnicity vs cross-ethnic pools.

The end product is a Word manuscript (`HLA_Registry_Size_CMIO_v*.docx`), regenerated each version. `.docx` files are intentionally **not** committed (see git log).

## Commands

```bash
pip install -r analysis/requirements.txt   # deps (build_report.py also needs python-docx)
bash analysis/run_all.sh                    # core pipeline: steps 01→05 + plot_coverage
pytest tests/ -v                            # 35 tests (run from repo root)
python build_report_v215c.py                # regenerate the manuscript .docx (run from repo root)
                                            # NOTE: builders are gitignored — local only
```

Individual analysis steps run from inside `analysis/` (e.g. `cd analysis && python 09_bootstrap_ci.py`). They read/write `analysis/data/*.csv` and `analysis/figures/*.png`.

## Architecture — the part that isn't obvious from filenames

**Two-tier pipeline.** `run_all.sh` only runs the *core* (01 ingest → 02 allele_freq → 03 hwe → 04 registry_model → plot_coverage → 05 report). Scripts **06–16 and the bootstrap are NOT in run_all.sh** — they are standalone, run manually, and produce the CSV/PNG inputs that `build_report.py` consumes. So `run_all.sh` passing ≠ the manuscript being current. To rebuild the paper end-to-end you must also run 06, 07, 09–16 before `build_report.py`.

**Numbered scripts are runners; logic lives in three library modules** (`hwe_test.py`, `registry_model.py`, `plot_coverage.py`) — these are what the tests import. Core model math in `registry_model.py`: `get_diplotype_frequencies` (HWE expansion), `compute_coverage` (C(N)=Σ f·[1−(1−f)^N]), `find_registry_size` (log-scale binary search), `get_combined_haplotype_freqs` (Singapore-weighted pool). EM haplotype phasing in `hwe_test.py:run_em_haplotypes` — **capped at 5,000 individuals/ethnicity** for speed (a known ~8% conservative bias, quantified by `15_em_convergence.py`).

**Data flows long-format.** `01_ingest.py` reads Excel/txt sources → `analysis/data/hla_clean.csv` (one row per sample×locus; `source` ∈ {`BMDP_OUT`,`SCBB_OUT`}). Downstream scripts filter on the `MAIN_SOURCES` constant. Ingestion is skipped if `hla_clean.csv` exists. Full input-format contract (column regexes, ethnicity codes, bypassing ingestion) is documented in `README.md` → "Using Your Own CMIO HLA Data".

**Two report generators at root — only one is current.** `build_report.py` is canonical: data-driven, reads `registry_size_ci.csv` + other CSVs, emits bootstrap-median registry sizes with 95% CIs. `generate_report.py` is an **older draft with hardcoded numbers** (different, smaller registry figures) — don't update it.

**The frequency floor sets the scale of every registry number.** v2.15c uses 1e-3; the live `analysis/data/` holds a 1e-6 re-run whose targets are ~2,098x larger and whose cross-group ranking is inverted. Comparative conclusions survive the change, absolute targets and rankings do not. Before quoting or editing any registry size, read `docs/explanation-frequency-floor.md`.

**Watch for two coexisting number sets.** `verification_summary.md`/`README` results tables and the `build_report.py` manuscript use different models (per-locus exact-match vs EM-phased bootstrap median), so e.g. "Chinese 95% 10/10" legitimately appears as both ~11.6k and ~42.8k. Don't "fix" one to match the other — confirm which model a table is reporting first.

## Conventions

- `VERSION.md` is the changelog and is updated per `.docx` release; bump it when regenerating the report.
- Figures are written by the script that owns the analysis (e.g. `10_ld_report.py` → `ld_heatmap_*.png`), then embedded by `build_report.py`.
- Reference PDFs and large `.xlsx` source data live untracked in the repo root.

---

# context-mode — MANDATORY routing rules

You have context-mode MCP tools available. These rules are NOT optional — they protect your context window from flooding. A single unrouted command can dump 56 KB into context and waste the entire session.

## BLOCKED commands — do NOT attempt these

### curl / wget — BLOCKED
Any Bash command containing `curl` or `wget` is intercepted and replaced with an error message. Do NOT retry.
Instead use:
- `ctx_fetch_and_index(url, source)` to fetch and index web pages
- `ctx_execute(language: "javascript", code: "const r = await fetch(...)")` to run HTTP calls in sandbox

### Inline HTTP — BLOCKED
Any Bash command containing `fetch('http`, `requests.get(`, `requests.post(`, `http.get(`, or `http.request(` is intercepted and replaced with an error message. Do NOT retry with Bash.
Instead use:
- `ctx_execute(language, code)` to run HTTP calls in sandbox — only stdout enters context

### WebFetch — BLOCKED
WebFetch calls are denied entirely. The URL is extracted and you are told to use `ctx_fetch_and_index` instead.
Instead use:
- `ctx_fetch_and_index(url, source)` then `ctx_search(queries)` to query the indexed content

## REDIRECTED tools — use sandbox equivalents

### Bash (>20 lines output)
Bash is ONLY for: `git`, `mkdir`, `rm`, `mv`, `cd`, `ls`, `npm install`, `pip install`, and other short-output commands.
For everything else, use:
- `ctx_batch_execute(commands, queries)` — run multiple commands + search in ONE call
- `ctx_execute(language: "shell", code: "...")` — run in sandbox, only stdout enters context

### Read (for analysis)
If you are reading a file to **Edit** it → Read is correct (Edit needs content in context).
If you are reading to **analyze, explore, or summarize** → use `ctx_execute_file(path, language, code)` instead. Only your printed summary enters context. The raw file content stays in the sandbox.

### Grep (large results)
Grep results can flood context. Use `ctx_execute(language: "shell", code: "grep ...")` to run searches in sandbox. Only your printed summary enters context.

## Tool selection hierarchy

1. **GATHER**: `ctx_batch_execute(commands, queries)` — Primary tool. Runs all commands, auto-indexes output, returns search results. ONE call replaces 30+ individual calls.
2. **FOLLOW-UP**: `ctx_search(queries: ["q1", "q2", ...])` — Query indexed content. Pass ALL questions as array in ONE call.
3. **PROCESSING**: `ctx_execute(language, code)` | `ctx_execute_file(path, language, code)` — Sandbox execution. Only stdout enters context.
4. **WEB**: `ctx_fetch_and_index(url, source)` then `ctx_search(queries)` — Fetch, chunk, index, query. Raw HTML never enters context.
5. **INDEX**: `ctx_index(content, source)` — Store content in FTS5 knowledge base for later search.

## Subagent routing

When spawning subagents (Agent/Task tool), the routing block is automatically injected into their prompt. Bash-type subagents are upgraded to general-purpose so they have access to MCP tools. You do NOT need to manually instruct subagents about context-mode.

## Output constraints

- Keep responses under 500 words.
- Write artifacts (code, configs, PRDs) to FILES — never return them as inline text. Return only: file path + 1-line description.
- When indexing content, use descriptive source labels so others can `ctx_search(source: "label")` later.

## ctx commands

| Command | Action |
|---------|--------|
| `ctx stats` | Call the `ctx_stats` MCP tool and display the full output verbatim |
| `ctx doctor` | Call the `ctx_doctor` MCP tool, run the returned shell command, display as checklist |
| `ctx upgrade` | Call the `ctx_upgrade` MCP tool, run the returned shell command, display as checklist |
