import argparse
import glob
import matplotlib.pyplot as plt
import numpy as np
import spinlip as sp


def main():
    parser = argparse.ArgumentParser(
        description="Equilibration analysis of the J1-J2 random-bond Ising model."
    )
    parser.add_argument("--dimension",              "-D", type=int,   default=2,                help="Lattice dimension (default: 2)")
    parser.add_argument("--linear_size",            "-L", type=int,   default=32,               help="Linear system size (default: 32)")
    parser.add_argument("--temperature",            "-T", type=float, default=0.908,            help="Temperature (default: 0.908)")
    parser.add_argument("--occupation_probability", "-P", type=float, default=0.5,             help="J2 bond occupation probability (default: 0.5)")
    parser.add_argument("--bond_strength",          "-J", type=float, default=0.1,             help="Weak coupling J2 (default: 0.1)")
    parser.add_argument("--disorder_seed",                type=int,   default=0,               help="Disorder seed for single-sample plot (default: 0)")
    parser.add_argument("--disorder_sample",              type=int,   default=0,               help="Disorder sample index for single-sample plot (default: 0)")
    parser.add_argument("--maximum_time",                 type=int,   default=100,             help="Number of sweeps to show (default: 250)")
    parser.add_argument("--path_in",                      type=str,   default="./../data/new/", help="Input data directory")
    parser.add_argument("--plot_path",                    type=str,   default="./../plots/",   help="Output directory for plots")
    args = parser.parse_args()


    DIMENSION      = args.dimension
    LINEAR_SIZE    = args.linear_size
    TEMPERATURE    = args.temperature
    P              = args.occupation_probability
    J2             = args.bond_strength
    DISORDERSEED   = args.disorder_seed
    DISORDERSAMPLE = args.disorder_sample
    MAXIMUM_TIME   = args.maximum_time
    PATH           = args.path_in
    PLOT_PATH      = args.plot_path

    NUM_SPINS = LINEAR_SIZE ** DIMENSION

    # --- single disorder sample ---
    pattern_rs = (PATH + f"J1J2Equilibration_D{DIMENSION}_L{LINEAR_SIZE}_P{P:.4f}_J{J2:.4f}_T{TEMPERATURE:.6f}"
                       + f"_DisorderSeed{DISORDERSEED}_DisorderSample{DISORDERSAMPLE}_*_RS_.dat")
    pattern_gs = (PATH + f"J1J2Equilibration_D{DIMENSION}_L{LINEAR_SIZE}_P{P:.4f}_J{J2:.4f}_T{TEMPERATURE:.6f}"
                       + f"_DisorderSeed{DISORDERSEED}_DisorderSample{DISORDERSAMPLE}_*_GS_.dat")
    file_names_rs = glob.glob(pattern_rs)
    file_names_gs = glob.glob(pattern_gs)

    data_rs          = np.loadtxt(file_names_rs[0])
    data_gs          = np.loadtxt(file_names_gs[0])
    time             = data_rs[:, 0]
    energy_rs        = data_rs[:, 1]
    magnetization_rs = np.abs(data_rs[:, 2])
    energy_gs        = data_gs[:, 1]
    magnetization_gs = np.abs(data_gs[:, 2])

    # --- all disorder samples ---
    pattern_all_rs = (PATH + f"J1J2Equilibration_D{DIMENSION}_L{LINEAR_SIZE}_P{P:.4f}_J{J2:.4f}_T{TEMPERATURE:.6f}"
                           + f"_*_RS_.dat")
    pattern_all_gs = (PATH + f"J1J2Equilibration_D{DIMENSION}_L{LINEAR_SIZE}_P{P:.4f}_J{J2:.4f}_T{TEMPERATURE:.6f}"
                           + f"_*_GS_.dat")
    file_names_all_rs = glob.glob(pattern_all_rs)
    file_names_all_gs = glob.glob(pattern_all_gs)
    print(f"number of RS files: {len(file_names_all_rs)}")
    print(f"number of GS files: {len(file_names_all_gs)}")

    NUM_SAMPLES_RS        = len(file_names_all_rs)
    NUM_SAMPLES_GS        = len(file_names_all_gs)
    data_energy_rs        = np.empty((NUM_SAMPLES_RS, MAXIMUM_TIME))
    data_magnetization_rs = np.empty((NUM_SAMPLES_RS, MAXIMUM_TIME))
    data_energy_gs        = np.empty((NUM_SAMPLES_GS, MAXIMUM_TIME))
    data_magnetization_gs = np.empty((NUM_SAMPLES_GS, MAXIMUM_TIME))

    for c, file_name in enumerate(file_names_all_rs):
        d = np.loadtxt(file_name)
        data_energy_rs[c, :]        = d[:MAXIMUM_TIME, 1] / NUM_SPINS
        data_magnetization_rs[c, :] = np.abs(d[:MAXIMUM_TIME, 2]) / NUM_SPINS
    for c, file_name in enumerate(file_names_all_gs):
        d = np.loadtxt(file_name)
        data_energy_gs[c, :]        = d[:MAXIMUM_TIME, 1] / NUM_SPINS
        data_magnetization_gs[c, :] = np.abs(d[:MAXIMUM_TIME, 2]) / NUM_SPINS

    avg_energy_rs     = np.empty(MAXIMUM_TIME)
    avg_energy_rs_err = np.empty(MAXIMUM_TIME)
    avg_energy_gs     = np.empty(MAXIMUM_TIME)
    avg_energy_gs_err = np.empty(MAXIMUM_TIME)
    avg_mag_rs        = np.empty(MAXIMUM_TIME)
    avg_mag_rs_err    = np.empty(MAXIMUM_TIME)
    avg_mag_gs        = np.empty(MAXIMUM_TIME)
    avg_mag_gs_err    = np.empty(MAXIMUM_TIME)

    for c in range(MAXIMUM_TIME):
        avg_energy_rs[c],  avg_energy_rs_err[c] = sp.calc_stats(data_energy_rs[:, c])
        avg_energy_gs[c],  avg_energy_gs_err[c] = sp.calc_stats(data_energy_gs[:, c])
        avg_mag_rs[c],     avg_mag_rs_err[c]    = sp.calc_stats(data_magnetization_rs[:, c])
        avg_mag_gs[c],     avg_mag_gs_err[c]    = sp.calc_stats(data_magnetization_gs[:, c])

    
    plt.close("all")
    plot_tag     = f"D{DIMENSION}_L{LINEAR_SIZE}_P{P:.4f}_J{J2:.4f}_T{TEMPERATURE:.6f}_DisorderSeed{DISORDERSEED}_DisorderSample{DISORDERSAMPLE}"
    plot_tag_avg = f"D{DIMENSION}_L{LINEAR_SIZE}_P{P:.4f}_J{J2:.4f}_T{TEMPERATURE:.6f}_DisorderAvg"
    t = time[:MAXIMUM_TIME]

        # single-sample figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.set_title(r'energy: $e$,  $L$ = '+str(LINEAR_SIZE)+', '+r'$T$ = '+str(TEMPERATURE))
    ax1.plot(time, energy_rs / NUM_SPINS, label="random start", color='red')
    ax1.plot(time, energy_gs / NUM_SPINS, label="ground state start", color='blue')
    ax1.set_xlabel("sweeps")
    ax1.set_ylabel(r"$e$")
    
    ax1.set_xlim([-10, MAXIMUM_TIME])
    ax1.legend()

    ax2.set_title(r'magnetization: $m$,  $L$ = '+str(LINEAR_SIZE)+', '+r'$T$ = '+str(TEMPERATURE))
    ax2.plot(time, magnetization_rs / NUM_SPINS, label="random start", color='red')
    ax2.plot(time, magnetization_gs / NUM_SPINS, label="ground state start", color='blue')
    ax2.set_xlabel("sweeps")
    ax2.set_ylabel(r"$|m|$")
    ax2.set_xlim([-10, MAXIMUM_TIME])
    ax2.legend()

    fig.tight_layout()
    fig.savefig(PLOT_PATH + f"J1J2Equilibration_{plot_tag}.pdf", bbox_inches="tight", pad_inches=0.03)

    # disorder-averaged figure
    fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(10, 4))
    ax3.set_title(r'disorder-avg energy: $[e]$,  $L$ = '+str(LINEAR_SIZE)+', '+r'$T$ = '+str(TEMPERATURE))
    ax3.plot(t, avg_energy_rs, label="random start", color='red')
    ax3.fill_between(t, avg_energy_rs - avg_energy_rs_err, avg_energy_rs + avg_energy_rs_err, alpha=0.3, color='red')
    ax3.plot(t, avg_energy_gs, label="ground state start", color='blue')
    ax3.fill_between(t, avg_energy_gs - avg_energy_gs_err, avg_energy_gs + avg_energy_gs_err, alpha=0.3, color='blue')
    ax3.set_xlabel("sweeps")
    ax3.set_ylabel(r"$e$")
    ax3.set_xlim([-10, MAXIMUM_TIME])
    ax3.legend()

    ax4.set_title(r'disorder-avg magnetization: $[|m|]$,  $L$ = '+str(LINEAR_SIZE)+', '+r'$T$ = '+str(TEMPERATURE))
    ax4.plot(t, avg_mag_rs, label="random start", color='red')
    ax4.fill_between(t, avg_mag_rs - avg_mag_rs_err, avg_mag_rs + avg_mag_rs_err, alpha=0.3, color='red')
    ax4.plot(t, avg_mag_gs, label="ground state start", color='blue')
    ax4.fill_between(t, avg_mag_gs - avg_mag_gs_err, avg_mag_gs + avg_mag_gs_err, alpha=0.3, color='blue')
    ax4.set_xlabel("sweeps")
    ax4.set_ylabel(r"$[|m|]$")
    ax4.set_xlim([-10, MAXIMUM_TIME])
    ax4.legend()

    fig2.tight_layout()
    fig2.savefig(PLOT_PATH + f"J1J2EquilibrationAvg_{plot_tag_avg}.pdf", bbox_inches="tight", pad_inches=0.03)
    plt.show()


if __name__ == "__main__":
    main()
