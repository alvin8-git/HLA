#!/bin/bash
set -e
cd "$(dirname "$0")"

echo "=== Step 1: Ingest data ==="
python3 01_ingest.py

echo "=== Step 2: Allele frequency verification ==="
python3 02_allele_freq.py

echo "=== Step 3: HWE tests + haplo.stats (Python EM) ==="
python3 03_hwe_test.py

echo "=== Step 4: Registry size model ==="
python3 04_registry_model.py

echo "=== Step 5: Coverage figures ==="
python3 plot_coverage.py

echo "=== Step 6: Report assembly ==="
python3 05_report.py

echo "Done. Outputs in analysis/data/ and analysis/figures/"
