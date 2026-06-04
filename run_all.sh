#!/bin/bash
set -e

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- 1. Clean output directories ---
echo "=== Cleaning data/new, data/results, data/bonds, plots ==="
rm -f "$ROOT_DIR/data/new/"*
rm -f "$ROOT_DIR/data/bonds/"*
rm -f "$ROOT_DIR/data/results/"*
rm -f "$ROOT_DIR/plots/"*

# --- 2. Compile ---
echo "=== Compiling ==="
cd "$ROOT_DIR/cplusplus"
make
cd "$ROOT_DIR"

# --- 3. Run simulation (fixed disorder, single sample) ---
echo "=== Running temperature scan ==="
cd "$ROOT_DIR/cplusplus"
bash run_temperature_scan_single_sample.sh
cd "$ROOT_DIR"

# --- 4. Disorder sampling analysis ---
echo "=== Running disorder sampling analysis ==="
cd "analysis"
MPLBACKEND=Agg python3 "j1j2_disorder_sampling_analysis.py" \
    --fixed_disorder True \
    --path_in "$ROOT_DIR/data/new/" \
    --path_out "$ROOT_DIR/data/results/" \
    --plot True \
    --save True

# --- 5. Comparison with exact solution ---
echo "=== Plotting comparison with exact solution ==="
MPLBACKEND=Agg python3 j1j2_comparison_exact_solution_plot.py
cd "$ROOT_DIR"

echo "=== Done. Results in data/results/, plots in data/results/ and plots/ ==="
