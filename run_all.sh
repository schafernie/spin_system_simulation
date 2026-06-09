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
cd "$ROOT_DIR/simulation"
make
cd "$ROOT_DIR"

# --- 3. Equilibration and autocorrelation ---
echo "=== Running equilibration ==="
cd "$ROOT_DIR/simulation"
bash run_equilibration.sh
cd "$ROOT_DIR"

echo "=== Running equilibration plot ==="
cd "$ROOT_DIR/scripts"
MPLBACKEND=Agg python3 j1j2_equilibration_plot.py \
    --path_in "$ROOT_DIR/data/new/" \
    --plot_path "$ROOT_DIR/plots/"
cd "$ROOT_DIR"

echo "=== Running autocorrelation plot ==="
cd "$ROOT_DIR/scripts"
MPLBACKEND=Agg python3 j1j2_autocorrelation_plot.py \
    --path_in "$ROOT_DIR/data/new/" \
    --plot_path "$ROOT_DIR/plots/"
cd "$ROOT_DIR"

# --- 4. Run simulation (fixed disorder, single sample) ---
echo "=== Running temperature scan ==="
cd "$ROOT_DIR/simulation"
bash run_temperature_scan_single_sample.sh
cd "$ROOT_DIR"

# --- 5. Disorder sampling analysis ---
echo "=== Running disorder sampling analysis ==="
cd "$ROOT_DIR/scripts"
MPLBACKEND=Agg python3 j1j2_disorder_sampling_analysis.py \
    --fixed_disorder True \
    --path_in "$ROOT_DIR/data/new/" \
    --path_out "$ROOT_DIR/data/results/" \
    --plot True \
    --save True
cd "$ROOT_DIR"

# --- 6. Comparison with exact solution ---
echo "=== Plotting comparison with exact solution ==="
cd "$ROOT_DIR/scripts"
MPLBACKEND=Agg python3 j1j2_comparison_exact_solution_plot.py
cd "$ROOT_DIR"

echo "=== Done. Results in data/results/, plots in data/results/ and plots/ ==="
