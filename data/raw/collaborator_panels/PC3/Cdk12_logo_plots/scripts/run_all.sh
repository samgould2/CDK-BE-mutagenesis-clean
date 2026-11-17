#!/usr/bin/env bash
# Runs the full pipeline end to end.
# Usage: bash run_all.sh   (run from inside the scripts/ directory)
set -e

echo "== Step 1: extract reference FASTA from annotated GenBank files =="
python3 01_extract_reference_fasta.py

echo
echo "== Step 2: align FASTQs and pileup codon-level base calls =="
python3 02_align_and_pileup.py

echo
echo "== Step 3: build amino-acid logo plots =="
python3 03_build_logo_plots.py

echo
echo "Done. See ../output/ for results."
