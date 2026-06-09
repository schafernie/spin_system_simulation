# Disordered Ising Model: Monte Carlo Simulation

Monte Carlo simulation of the J1-J2 random-bond Ising model on a $d$-dimensional hypercubic lattice. Each bond is independently drawn as a weak coupling $0<J_2<1$ with probability $p$ or a strong coupling $J_1 = 1$ with probability $1-p$. The Hamiltonian is

$$\mathcal{H} = -\sum_{\langle i,j \rangle} J_{ij} s_i s_j, \quad s_i = \pm 1,$$

where the quenched disorder $\{J_{ij}\}$ is sampled once per disorder realization and then held fixed during thermal sampling. The simulation uses the **Swendsen-Wang cluster algorithm** for efficient sampling near the critical temperature.

## Overview

The workflow has three stages:

1. **Build** — compile the C++ simulator with `make`
2. **Run** — generate disorder-sampling data with `mc_sim`
3. **Analyze** — compute disorder-averaged observables and Binder cumulant with the Python scripts in `scripts/`

## Quickstart

```bash
# Build
cd simulation
make

# Run everything (equilibration, temperature scan, analysis)
bash run_all.sh

# Or run steps individually:
bash run_equilibration.sh                        # check equilibration (produces plots in plots/)
bash run_temperature_scan_single_sample.sh       # temperature scan
cd ../scripts
python3 j1j2_disorder_sampling_analysis.py       # disorder-averaged observables
```

Default parameters: lattice size $L=32$, dimension $d=2$, bond occupation probability $p=0.5$, weak coupling $J_2=0.1$.

## File Structure

```
spin_system_simulation/
├── README.md
├── LICENSE
├── requirements.txt
├── simulation/
│   ├── Makefile                               # Build system (supports GCC and Intel icpx/MKL)
│   ├── testrun.sh                             # Quick test: single sampling run
│   ├── run_equilibration.sh                   # Equilibration run across disorder samples
│   ├── run_temperature_scan.sh                # Temperature scan, disorder seed increments with T
│   ├── run_temperature_scan_single_sample.sh  # Temperature scan, fixed disorder seed
│   └── src/
│       ├── appl/main.cpp                      # Entry point and CLI argument parsing
│       ├── simulation/                        # Top-level simulation routines
│       ├── mcmc/                              # Swendsen-Wang MCMC
│       ├── interaction/                       # J1-J2 bond disorder
│       ├── ising_configuration/               # Spin configuration
│       ├── lattice/                           # Hypercubic lattice geometry
│       ├── clusters/                          # Cluster data structure for Swendsen-Wang
│       ├── observables/                       # Energy and magnetization
│       ├── random_numbers/                    # RNG wrapper
│       └── statistics/                        # Moment estimators
├── scripts/
│   ├── j1j2_equilibration_plot.py            # Equilibration convergence plots (single sample + disorder average)
│   ├── j1j2_disorder_sampling_analysis.py    # Disorder-averaged magnetization, energy, Binder cumulant
│   └── j1j2_comparison_exact_solution_plot.py # Compare MC results against exact free energy
└── data/                                      # Generated simulation data
```

## Building

```bash
cd simulation
make           # default: GCC + OpenMP
make USE_MKL=1 # Intel icpx + MKL parallel RNG
make clean     # remove objects and binary
```

The binary `mc_sim` is placed in `simulation/`.

## Running the Simulator

```bash
./mc_sim -O <mode> -L <L> -D <D> -T <T> -P <p> -J <J2> [options]
```

### Modes

| `-O` | Mode | Description |
|---|---|---|
| `0` | Equilibration | Records energy and magnetization sweep-by-sweep from two starting conditions (random and all-up ground state) |
| `1` | Disorder sampling | Samples thermal observables for multiple disorder realizations and writes per-realization moments |
| `-1` | Debugging | Prints lattice, RNG, bond, and MCMC diagnostics |

### Parameters

| Argument | Default | Description |
|---|---|---|
| `-L <int>` | 4 | Linear system size |
| `-D <int>` | 2 | Lattice dimension |
| `-T <double>` | 1.0 | Temperature |
| `-P <double>` | 0.5 | Probability of a weak $J_2$ bond |
| `-J <double>` | 0.1 | Weak coupling strength $J_2$ (strong coupling $J_1 = 1$) |
| `--DisorderSeed <int>` | 0 | Seed for the disorder RNG |
| `--DisorderSample <int>` | 0 | Starting disorder-realization index |
| `--NrDisorderSamples <int>` | 1 | Number of disorder realizations to sample |
| `--EquilibrationTime <int>` | 1000 | Swendsen-Wang sweeps before recording |
| `--ThermalSeed <int>` | 0 | Seed for the thermal RNG |
| `--ThermalSamples <int>` | 100 | Thermal samples per disorder realization |
| `--TimeBetweenSamples <int>` | 10 | Swendsen-Wang sweeps between consecutive samples |
| `--FixedDisorder` | off | Keep the same bond configuration across all disorder iterations |

### Example commands

```bash
# Disorder sampling at a single temperature
./mc_sim -O 1 -L 32 -D 2 -T 0.90 -J 0.1 -P 0.5 \
    --DisorderSeed 0 --DisorderSample 0 --NrDisorderSamples 10 \
    --EquilibrationTime 1000 --ThermalSeed 0 --ThermalSamples 1000 --TimeBetweenSamples 25

# Temperature scan across T = 0.70 … 1.20 in steps of 0.05
bash run_temperature_scan.sh

# Fixed-disorder temperature scan (same bond configuration at every T)
bash run_temperature_scan_single_sample.sh
```

## Output Files

All output is written to `data/new/` (created automatically if absent).

**Disorder sampling** (`-O 1`):

- `J1J2DisorderSampling_D{D}_L{L}_P{p}_J{J2}_T{T}_ThermalSamples{n}_DisorderSeed{s}_DisorderSample{c}_.dat`
- `J1J2DisorderSamplingFixed_D{D}_L{L}_P{p}_J{J2}_T{T}_ThermalSamples{n}_DisorderSeed{s}_DisorderSample{c}_.dat` (with `--FixedDisorder`)

Each row is one disorder realization with columns:
`disorder_sample  |m|  |m|^2  |m|^3 |m|^4  e  e^2  e^3  e^4`

where $m = M/N$ and $e = E/N$ are per-spin quantities.

## Analysis Scripts

All scripts read from `../data/new/` and write results to `../data/results/`.

**`scripts/j1j2_equilibration_plot.py`** — equilibration convergence

Plots energy and magnetization sweep-by-sweep for a single disorder sample and the disorder average (with error bands) from two starting conditions (random and ground state), to verify that the system has equilibrated before sampling.

**`scripts/j1j2_disorder_sampling_analysis.py`** — disorder-averaged thermodynamics

Loads disorder-sampling output, computes disorder averages of $\langle|m|\rangle$, $\langle m^2\rangle$, $\langle m^4\rangle$, $\langle e\rangle$, and the Binder cumulant $U_4 = 1 - \langle m^4\rangle / (3\langle m^2\rangle^2)$ with bootstrap error bars, and saves results and plots as a function of temperature.

Key options: `--fixed_disorder` to load fixed-disorder files, `--disorder_seed` / `--disorder_sample` to select the realization, `--plot` / `--save` to control output.

**`scripts/j1j2_autocorrelation_plot.py`** — autocorrelation analysis

Computes the normalized autocorrelation function of the magnetization and energy from equilibration runs, averaging over disorder samples, and plots both the direct estimator and the FFT estimator alongside an exponential fit to extract the autocorrelation time $\tau_\mathrm{exp}$.

**`scripts/j1j2_comparison_exact_solution_plot.py`** — comparison with exact solution 


Loads the exact free energy from `data/exact/` and compares energy and specific heat (derived via numerical differentiation) against the MC disorder averages. The numerical exact solution is taken from Phys. Rev. E, 87, 043303 (arXiv:1301.1252v1), see https://github.com/a-alan-middleton/IsingPartitionFn.


## Dependencies

**C++**: C++17 compiler, OpenMP. Optionally Intel oneAPI (`icpx`, MKL) for parallel RNG (`USE_MKL=1`).

**Python**: NumPy, Matplotlib. Install with:

```bash
pip install -r requirements.txt
```

## License

MIT License — see [LICENSE](LICENSE) for details.

Copyright (C) 2026 Lambert Münster
