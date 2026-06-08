import numpy as np
import matplotlib.pyplot as plt

# --- Parameters ---
DATA_PATH = "./../data/"
PLOT_PATH = "./../plots/"

LINEAR_SIZE = 32
P = 0.5
J = 0.1
DISORDER_SEED = 0
DISORDER_SAMPLE = 0

NUM_SITES = LINEAR_SIZE ** 2

# --- Load exact free energy ---
exact_file = (
    f"J1J2FreeEnergy_P{P:.4f}_J{J:.4f}"
    f"_DisorderSeed{DISORDER_SEED}_DisorderSample{DISORDER_SAMPLE}_MKL_.dat"
)
exact_data = np.loadtxt(DATA_PATH + "exact/" + exact_file, skiprows=1)
T_exact = exact_data[:, 0]
f = exact_data[:, 1] / NUM_SITES 

# Thermodynamic quantities from free energy via numerical differentiation
df_dT = np.gradient(f, T_exact)
e_exact = f - T_exact * df_dT

df2_dT2 = np.gradient(df_dT, T_exact)
mask = slice(5, -5)  # trim edges where gradient is noisy
T_c = T_exact[mask]
c_exact = -T_c * df2_dT2[mask]

# --- Load MC energy ---
mc_energy_file = f"J1J2E_D2_L{LINEAR_SIZE}_P{P:.4f}_J{J:.4f}_DisorderSeed{DISORDER_SEED}_DisorderSample{DISORDER_SAMPLE}_.dat"
mc_energy_data = np.loadtxt(DATA_PATH + "results/" + mc_energy_file, skiprows=1)
T_mc_e = mc_energy_data[:, 0]
e_mc = mc_energy_data[:, 1]
e_mc_err = mc_energy_data[:, 2]



#%%
# --- Load MC energy variance (for specific heat) ---
mc_var_file = f"J1J2VarE_D2_L{LINEAR_SIZE}_P{P:.4f}_J{J:.4f}_DisorderSeed{DISORDER_SEED}_DisorderSample{DISORDER_SAMPLE}_.dat"
mc_var_data = np.loadtxt(DATA_PATH + "results/" + mc_var_file, skiprows=1)
T_mc_c = mc_var_data[:, 0]
c_mc = mc_var_data[:, 1] * NUM_SITES  / T_mc_c ** 2
c_mc_err = mc_var_data[:, 2] * NUM_SITES  / T_mc_c ** 2

# --- Plot ---
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
fig.suptitle(
    rf"J1-J2 disordered Ising model: $L={LINEAR_SIZE}$, $J={J}$, "
    rf"$P={P}$"
)

# Free energy
axes[0].plot(T_exact, f, ls="--", label="exact")
axes[0].set_xlabel("T")
axes[0].set_ylabel("free energy per site")
axes[0].legend()

# Energy
axes[1].plot(T_exact, e_exact, ls="--", label="exact")
axes[1].errorbar(T_mc_e, e_mc, e_mc_err, marker="x", color="red", ls=" ", label="MC")
axes[1].set_xlabel("T")
axes[1].set_ylabel("energy per site")
axes[1].legend()

# Specific heat
axes[2].plot(T_c, c_exact, ls="--", label="exact")
axes[2].errorbar(T_mc_c, c_mc, c_mc_err, marker="x", color="red", ls=" ", label="MC")
axes[2].set_xlabel("T")
axes[2].set_ylabel("specific heat per site")
axes[2].legend()

fig.tight_layout()

save_name = f"J1J2Comparison_L{LINEAR_SIZE}_P{P:.4f}_J{J:.4f}_DisorderSeed{DISORDER_SEED}_DisorderSample{DISORDER_SAMPLE}.pdf"
fig.savefig(PLOT_PATH + save_name, bbox_inches="tight", pad_inches=0.03)
plt.show()
