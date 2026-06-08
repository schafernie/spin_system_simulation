import argparse
import glob
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import spinlip as sp



def main():
    parser = argparse.ArgumentParser(
        description="Disorder-averaged thermodynamics of the J1-J2 random-bond Ising model."
    )
    parser.add_argument("--dimension",       "-D", type=int,   default=2,    help="Lattice dimension (default: 2)")
    parser.add_argument("--linear_size",     "-L", type=int,   default=32,   help="Linear system size (default: 32)")
    parser.add_argument("--thermal_samples",       type=int,   default=1000, help="Thermal samples per disorder realization (default: 1000)")
    parser.add_argument("--occupation_probability", "-P", type=float, default=0.5,  help="J2 bond occupation probability (default: 0.5)")
    parser.add_argument("--bond_strength",       "-J", type=float, default=0.1,  help="Weak coupling J2 (default: 0.1)")
    parser.add_argument("--nr_boot_samples",       type=int,   default=250,  help="Bootstrap samples for Binder cumulant (default: 250)")
    parser.add_argument("--rng_seed",              type=int,   default=None, help="RNG seed for bootstrap (default: linear_size)")
    parser.add_argument("--path_in",               type=str,   default="./../data/new/", help="Input data directory")
    parser.add_argument("--path_out",              type=str,   default="./../data/results/",  help="Output directory (used with --save)")
    parser.add_argument("--plot_path",             type=str,   default="./../plots/",         help="Output directory for plots")
    parser.add_argument("--fixed_disorder",          type=bool,  default=True,   help="Load fixed-disorder files (J1J2DisorderSamplingFixed)")
    parser.add_argument("--disorder_seed",           type=int,   default=0,    help="Disorder seed (used with --fixed_disorder)")
    parser.add_argument("--disorder_sample",         type=int,   default=0,    help="Disorder sample (used with --fixed_disorder)")
    parser.add_argument("--save",                    type=bool,  default=True,    help="Save results to disk")
    parser.add_argument("--plot",                    type=bool,  default=True,    help="Save plots")
    args = parser.parse_args()

    DIMENSION       = args.dimension
    LINEAR_SIZE     = args.linear_size
    THERMAL_SAMPLES = args.thermal_samples
    P2              = args.occupation_probability
    J2              = args.bond_strength
    NR_BOOT_SAMPLES = args.nr_boot_samples
    RNG_SEED        = args.rng_seed if args.rng_seed is not None else LINEAR_SIZE
    PATH_IN         = args.path_in
    PATH_OUT        = args.path_out
    PLOT_PATH       = args.plot_path
    DO_PLOT         = args.plot
    DO_SAVE         = args.save
    FIXED_DISORDER  = args.fixed_disorder
    DISORDER_SEED   = args.disorder_seed
    DISORDER_SAMPLE = args.disorder_sample

    NR_SITES        = LINEAR_SIZE ** DIMENSION
    file_prefix     = "J1J2DisorderSamplingFixed" if FIXED_DISORDER else "J1J2DisorderSampling"
    disorder_tag    = f"_DisorderSeed{DISORDER_SEED}_DisorderSample{DISORDER_SAMPLE}" if FIXED_DISORDER else ""
    
    rng = np.random.default_rng(RNG_SEED)

    # --- discover temperatures from available files ---
    pattern = (
        f"{file_prefix}_D{DIMENSION}_L{LINEAR_SIZE}"
        f"_P{P2:.4f}_J{J2:.4f}"
        f"_*_ThermalSamples{THERMAL_SAMPLES}_*_.dat"
    )
    all_files = glob.glob(PATH_IN + pattern)
    print(f"Pattern: {PATH_IN + pattern}")
    print(f"Files found: {len(all_files)}")

    temperatures = []
    for fn in all_files:
        with open(fn) as f:
            line = f.readline()
        match = re.search(r"(?<=temperature: )([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)", line)
        temperatures.append(float(match[0]))
    temperatures = np.array(sorted(set(temperatures)))
    print(f"Temperatures: {temperatures}")

    num_temperatures   = len(temperatures)
    m          = np.empty(num_temperatures)
    m_err      = np.empty(num_temperatures)
    m2         = np.empty(num_temperatures)
    m2_err     = np.empty(num_temperatures)
    m4         = np.empty(num_temperatures)
    m4_err     = np.empty(num_temperatures)
    e          = np.empty(num_temperatures)
    e_err      = np.empty(num_temperatures)
    m_var      = np.empty(num_temperatures)
    m_var_err  = np.empty(num_temperatures)
    e_var      = np.empty(num_temperatures)
    e_var_err  = np.empty(num_temperatures)
    binder     = np.empty(num_temperatures)
    binder_err = np.empty(num_temperatures)
    nsamp      = np.empty(num_temperatures)
    boot       = np.empty(NR_BOOT_SAMPLES)

    # --- load data and compute observables ---
    for t, temp in enumerate(temperatures):
        pattern_t = (
            f"{file_prefix}_D{DIMENSION}_L{LINEAR_SIZE}"
            f"_P{P2:.4f}_J{J2:.4f}"
            f"_T{temp:.6f}*_ThermalSamples{THERMAL_SAMPLES}_*_.dat"
        )
        file_names = glob.glob(PATH_IN + pattern_t)
        print(f"T={temp:.6f}: {len(file_names)} file(s)")

        # columns: disorder_sample:  |m| ,  m^2 , m^3 , m^4 , e , e^2 , e^3 , e^4
        data  = np.loadtxt(file_names[0], skiprows=2)
        mall  = data[:, 1]
        m2all = data[:, 2]
        m4all = data[:, 4]
        eall  = data[:, 5]
        e2all = data[:, 6]
        for fn in file_names[1:]:
            data  = np.loadtxt(fn, skiprows=2)
            mall  = np.append(mall,  data[:, 1])
            m2all = np.append(m2all, data[:, 2])
            m4all = np.append(m4all, data[:, 4])
            eall  = np.append(eall,  data[:, 5])
            e2all = np.append(e2all, data[:, 6])

        nsamp[t] = len(mall)
        print(f"  disorder samples: {int(nsamp[t])}")

        m[t],     m_err[t]     = sp.calc_stats(mall)
        m2[t],    m2_err[t]    = sp.calc_stats(m2all)
        m4[t],    m4_err[t]    = sp.calc_stats(m4all)
        e[t],     e_err[t]     = sp.calc_stats(eall)
        m_var[t], m_var_err[t] = sp.calc_stats(
            THERMAL_SAMPLES / (THERMAL_SAMPLES - 1) * (m2all - mall**2)
        )
        e_var[t], e_var_err[t] = sp.calc_stats(
            THERMAL_SAMPLES / (THERMAL_SAMPLES - 1) * (e2all - eall**2)
        )
        # Binder cumulant B = 1 - <m^4> / (3 * <m^2>^2)
        binder[t] = 1 - np.mean(m4all) / (3 * np.mean(m2all)**2)

        n   = int(nsamp[t])
        idx = np.arange(n)
        for b in range(NR_BOOT_SAMPLES):
            ri      = rng.choice(idx, n, replace=True)
            boot[b] = 1 - np.mean(m4all[ri]) / (3 * np.mean(m2all[ri])**2)
        binder[t]     = 2 * binder[t] - np.mean(boot)
        binder_err[t] = np.std(boot)

    # --- save ---
    if DO_SAVE:
        os.makedirs(PATH_OUT, exist_ok=True)
        tag = f"D{DIMENSION}_L{LINEAR_SIZE}_P{P2:.4f}_J{J2:.4f}{disorder_tag}"
        np.savetxt(PATH_OUT + f"J1J2MagMoments_{tag}_.dat",
                   np.vstack((temperatures, m, m_err, m2, m2_err, m4, m4_err, nsamp)).T,
                   header="temperatures , |m| , err , m2 , err , m4 , err , nr. samples")
        np.savetxt(PATH_OUT + f"J1J2VarM_{tag}_.dat",
                   np.vstack((temperatures, m_var, m_var_err, nsamp)).T,
                   header="temperatures , var(|m|) , err , nr. samples")
        np.savetxt(PATH_OUT + f"J1J2VarE_{tag}_.dat",
                   np.vstack((temperatures, e_var, e_var_err, nsamp)).T,
                   header="temperatures , var(e) , err , nr. samples")
        np.savetxt(PATH_OUT + f"J1J2E_{tag}_.dat",
                   np.vstack((temperatures, e, e_err, nsamp)).T,
                   header="temperatures , e , err , nr. samples")
        np.savetxt(PATH_OUT + f"J1J2Binder_{tag}_.dat",
                   np.vstack((temperatures, binder, binder_err, nsamp)).T,
                   header="temperatures , U4 , err , nr. samples  (U4 = 1 - <m^4>/(3*<m^2>^2))")

    # --- plot ---
    if DO_PLOT:
        plt.close("all")
        plot_tag = f"P{P2:.4f}_J{J2:.4f}{disorder_tag}"

        fig, ax = plt.subplots()
        ax.set_title("energy per spin")
        ax.errorbar(temperatures, e, e_err, ls=" ", marker="x")
        ax.set_xlabel(r"$T$")
        ax.set_ylabel(r"$e$")
        fig.savefig(PLOT_PATH + f"J1J2E_{plot_tag}.pdf", bbox_inches="tight", pad_inches=0.03)

        fig, ax = plt.subplots()
        ax.set_title("magnetization per spin")
        ax.errorbar(temperatures, m, m_err, ls=" ", marker="x")
        ax.set_xlabel(r"$T$")
        ax.set_ylabel(r"$|m|$")
        fig.savefig(PLOT_PATH + f"J1J2M_{plot_tag}.pdf", bbox_inches="tight", pad_inches=0.03)

        fig, ax = plt.subplots()
        ax.set_title("variance of magnetization per spin")
        ax.errorbar(temperatures, m_var, m_var_err, ls=" ", marker="x")
        ax.set_xlabel(r"$T$")
        ax.set_ylabel(r"$\mathrm{var}(|m|)$")
        fig.savefig(PLOT_PATH + f"J1J2VarM_{plot_tag}.pdf", bbox_inches="tight", pad_inches=0.03)

        fig, ax = plt.subplots()
        ax.set_title("variance of energy per spin")
        ax.errorbar(temperatures, e_var, e_var_err, ls=" ", marker="x")
        ax.set_xlabel(r"$T$")
        ax.set_ylabel(r"$\mathrm{var}(e)$")
        fig.savefig(PLOT_PATH + f"J1J2VarE_{plot_tag}.pdf", bbox_inches="tight", pad_inches=0.03)

        fig, ax = plt.subplots()
        ax.set_title("square magnetization per spin")
        ax.errorbar(temperatures, NR_SITES * m2, NR_SITES * m2_err, ls=" ", marker="o")
        ax.set_xlabel(r"$T$")
        ax.set_ylabel(r"$\langle M^2 \rangle / N$")
        fig.savefig(PLOT_PATH + f"J1J2M2_{plot_tag}.pdf", bbox_inches="tight", pad_inches=0.03)

        fig, ax = plt.subplots()
        ax.set_title("Binder cumulant")
        ax.errorbar(temperatures, binder, binder_err, ls=" ", marker="o")
        ax.set_xlabel(r"$T$")
        ax.set_ylabel(r"$U_4 = 1 - \langle m^4 \rangle / (3\langle m^2 \rangle^2)$")
        fig.savefig(PLOT_PATH + f"J1J2Binder_{plot_tag}.pdf", bbox_inches="tight", pad_inches=0.03)

        plt.show()


if __name__ == "__main__":
    main()
