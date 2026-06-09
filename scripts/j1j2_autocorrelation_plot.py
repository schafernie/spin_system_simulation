import glob
import argparse
import matplotlib.pyplot as plt
import numpy as np
import spinlip as sp


def autocorrelation(q, t_max):
    autocorr = np.zeros(t_max)
    for t in range(t_max):
        autocorr[t] = chi(q, t)
    return autocorr


def chi(x, delta):
    n = len(x)
    return (np.mean(x[0:(n - delta)] * np.roll(x, -delta, axis=0)[0:(n - delta)])
            - np.mean(x[0:(n - delta)]) * np.mean(x[delta:]))

def main():
    parser = argparse.ArgumentParser(
        description="Equilibration analysis of the J1-J2 random-bond Ising model."
    )
    parser.add_argument("--dimension",              "-D", type=int,   default=2,                help="Lattice dimension (default: 2)")
    parser.add_argument("--linear_size",            "-L", type=int,   default=32,               help="Linear system size (default: 32)")
    parser.add_argument("--temperature",            "-T", type=float, default=0.905883,           help="Temperature (default: 0.905883)")
    parser.add_argument("--occupation_probability", "-P", type=float, default=0.5,              help="J2 bond occupation probability (default: 0.5)")
    parser.add_argument("--bond_strength",          "-J", type=float, default=0.1,              help="Weak coupling J2 (default: 0.1)")
    parser.add_argument("--t_eq",                         type=int,   default=1000,             help="Equilibration time to skip (default: 1000)")
    parser.add_argument("--t_max",                        type=int,   default=50,               help="Maximum autocorrelation lag (default: 50)")
    parser.add_argument("--t_int",                        type=int,   default=25,               help="Integration cutoff for tau_int (default: 25)")
    parser.add_argument("--path_in",                      type=str,   default="./../data/new/", help="Input data directory")
    parser.add_argument("--plot_path",                    type=str,   default="./../plots/",    help="Output directory for plots")
    args = parser.parse_args([])

    DIMENSION   = args.dimension
    LINEAR_SIZE = args.linear_size
    TEMPERATURE = args.temperature
    P           = args.occupation_probability
    J2          = args.bond_strength
    T_EQ        = args.t_eq
    T_MAX       = args.t_max
    T_INT       = args.t_int
    PATH_IN     = args.path_in
    PLOT_PATH   = args.plot_path

    NUM_SPINS = LINEAR_SIZE ** DIMENSION
    COL_E = 1  # energy column
    COL_M = 2  # magnetization column
    FIT_MIN = 2
    FIT_MAX = 10

    pattern = (PATH_IN + f"J1J2Equilibration_D{DIMENSION}_L{LINEAR_SIZE}_P{P:.4f}_J{J2:.4f}_T{TEMPERATURE:.6f}"
                       + "_*_RS_.dat")
    file_names = glob.glob(pattern)

    print("Autocorrelation: ")
    print("First Data File: " + file_names[0])

    n_files = len(file_names)

    autocorr_samples_m     = np.zeros([T_MAX, n_files])
    autocorr_fft_samples_m = np.zeros([T_MAX, n_files])
    tau_int_samples_m      = np.zeros(n_files)

    autocorr_samples_e     = np.zeros([T_MAX, n_files])
    autocorr_fft_samples_e = np.zeros([T_MAX, n_files])
    tau_int_samples_e      = np.zeros(n_files)

    for c in range(n_files):
        data = np.loadtxt(file_names[c], skiprows=3)

        # magnetization
        magnetization = np.abs(data[T_EQ:, COL_M]) / NUM_SPINS
        autocorr_samples_m[:, c] = autocorrelation(magnetization, T_MAX)
        tau_int_samples_m[c] = np.sum(autocorr_samples_m[1:T_INT, c] / autocorr_samples_m[0, c]) + 0.5
        magnetization = magnetization - np.mean(magnetization)
        magnetization = np.abs(np.fft.fft(magnetization)) ** 2
        autocorr_fft_samples_m[:, c] = np.abs(np.fft.ifft(magnetization))[:T_MAX]

        # energy
        energy = data[T_EQ:, COL_E] / NUM_SPINS
        autocorr_samples_e[:, c] = autocorrelation(energy, T_MAX)
        tau_int_samples_e[c] = np.sum(autocorr_samples_e[1:T_INT, c] / autocorr_samples_e[0, c]) + 0.5
        energy = energy - np.mean(energy)
        energy = np.abs(np.fft.fft(energy)) ** 2
        autocorr_fft_samples_e[:, c] = np.abs(np.fft.ifft(energy))[:T_MAX]

    autocorr_m         = np.zeros(T_MAX)
    autocorr_err_m     = np.zeros(T_MAX)
    autocorr_fft_m     = np.zeros(T_MAX)
    autocorr_fft_err_m = np.zeros(T_MAX)

    autocorr_e         = np.zeros(T_MAX)
    autocorr_err_e     = np.zeros(T_MAX)
    autocorr_fft_e     = np.zeros(T_MAX)
    autocorr_fft_err_e = np.zeros(T_MAX)

    for c in range(T_MAX):
        autocorr_m[c],     autocorr_err_m[c]     = sp.calc_stats(autocorr_samples_m[c, :])
        autocorr_fft_m[c], autocorr_fft_err_m[c] = sp.calc_stats(autocorr_fft_samples_m[c, :])
        autocorr_e[c],     autocorr_err_e[c]     = sp.calc_stats(autocorr_samples_e[c, :])
        autocorr_fft_e[c], autocorr_fft_err_e[c] = sp.calc_stats(autocorr_fft_samples_e[c, :])

    autocorr_err_m     = autocorr_err_m     / autocorr_m[0]
    autocorr_m         = autocorr_m         / autocorr_m[0]
    autocorr_fft_err_m = autocorr_fft_err_m / autocorr_fft_m[0]
    autocorr_fft_m     = autocorr_fft_m     / autocorr_fft_m[0]

    autocorr_err_e     = autocorr_err_e     / autocorr_e[0]
    autocorr_e         = autocorr_e         / autocorr_e[0]
    autocorr_fft_err_e = autocorr_fft_err_e / autocorr_fft_e[0]
    autocorr_fft_e     = autocorr_fft_e     / autocorr_fft_e[0]

    tau_int_m, tau_int_err_m = sp.calc_stats(tau_int_samples_m)
    tau_int_e, tau_int_err_e = sp.calc_stats(tau_int_samples_e)

  

    sweeps = np.arange(0, T_MAX)

    print(f"Integrated autocorrelation time (m): {tau_int_m:.4f} +/- {tau_int_err_m:.4f}")
    print(f"Integrated autocorrelation time (e): {tau_int_e:.4f} +/- {tau_int_err_e:.4f}")
    fname_suffix = f"D{DIMENSION}_L{LINEAR_SIZE}_P{P:.4f}_J{J2:.4f}_T{TEMPERATURE:.6f}.pdf"

    def fit_func(r, a, b):
        return a * np.exp(-r / b)

    plt.close("all")

    _, (ax_m, ax_e) = plt.subplots(1, 2, figsize=(12, 5))

    ax_m.set_title(r"Autocorrelation: magnetization")
    ax_m.set_yscale('log')
    ax_m.errorbar(sweeps, autocorr_m, autocorr_err_m, marker='x', color="green", ls='--', label=r'sum estimator', zorder=1)
    ax_m.plot(sweeps, autocorr_fft_m, color="black", label=r'FFT estimator')
    ax_m.axhline(y=1/np.e, color="black", ls='--', label=r'$1/e$')
    try:
        p_m, p_err_m, _ = sp.fitting_routine(fit_func, sweeps[FIT_MIN:FIT_MAX], autocorr_m[FIT_MIN:FIT_MAX],
                                             [1., 1.], autocorr_err_m[FIT_MIN:FIT_MAX])
        ax_m.plot(sweeps[FIT_MIN:FIT_MAX], fit_func(sweeps[FIT_MIN:FIT_MAX], *p_m), color="red", lw=2, zorder=2,
                  label=rf'fit $\tau_{{exp}}={p_m[1]:.1f} \pm {p_err_m[1]:.1f}$')
        print(f"Exponential autocorrelation time (m): {p_m[1]:.4f} +/- {p_err_m[1]:.4f}")
    except RuntimeError:
        print("Fit failed for magnetization autocorrelation.")
    ax_m.set_xlabel("mc sweep")
    ax_m.set_ylabel(r"$c$")
    ax_m.legend()

    ax_e.set_title(r"Autocorrelation: energy")
    ax_e.set_yscale('log')
    ax_e.errorbar(sweeps, autocorr_e, autocorr_err_e, marker='x', color="blue", ls='--', label=r'sum estimator', zorder=1)
    ax_e.plot(sweeps, autocorr_fft_e, color="black", label=r'FFT estimator')
    ax_e.axhline(y=1/np.e, color="black", ls='--', label=r'$1/e$')
    try:
        p_e, p_err_e, _ = sp.fitting_routine(fit_func, sweeps[FIT_MIN:FIT_MAX], autocorr_e[FIT_MIN:FIT_MAX],
                                             [1., 1.], autocorr_err_e[FIT_MIN:FIT_MAX])
        ax_e.plot(sweeps[FIT_MIN:FIT_MAX], fit_func(sweeps[FIT_MIN:FIT_MAX], *p_e), color="red", lw=2, zorder=2,
                  label=rf'fit $\tau_{{exp}}={p_e[1]:.1f} \pm {p_err_e[1]:.1f}$')
        print(f"Exponential autocorrelation time (e): {p_e[1]:.4f} +/- {p_err_e[1]:.4f}")
    except RuntimeError:
        print("Fit failed for energy autocorrelation.")
    ax_e.set_xlabel("mc sweep")
    ax_e.set_ylabel(r"$c$")
    ax_e.legend()

    plt.tight_layout()
    plt.savefig(PLOT_PATH + "Autocorrelation_" + fname_suffix)
    plt.show()


if __name__ == "__main__":
    main()
