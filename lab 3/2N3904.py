# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.axisbelow'] = True
%config InlineBackend.figure_format = 'retina'
# %%
small_sig = pd.read_csv("data/small_sig.csv")
vcesat = pd.read_csv("data/vcesat.csv")
# %%
phi_th = 0.026 # [V]
# %%
plt.grid()
plt.scatter(small_sig["I_C"], small_sig["I_C"]/phi_th, label="calculated")
plt.scatter(small_sig["I_C"], small_sig["ic"]/small_sig["vb"], label="measured")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("$I_C$")
plt.ylabel("$g_m$")
plt.legend()
plt.title("$g_m$ vs $I_C$")
# %%
plt.grid()
plt.scatter(small_sig["I_C"], small_sig["I_C"]/small_sig["I_B"])
plt.xscale("log")
plt.xlabel("$I_C$")
plt.ylabel("$\\beta$")
plt.ylim(0, 400)
plt.title("$\\beta$ vs $I_C$")
# %%
plt.grid()
plt.scatter(small_sig["I_C"], phi_th/small_sig["I_B"], label="calculated")
plt.scatter(small_sig["I_C"], small_sig["vb"]/small_sig["ib"], label="measured")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("$I_C$")
plt.ylabel("$r_\\pi$")
plt.legend()
plt.title("$r_\\pi$ vs $I_C$")
# %%
plt.grid()
plt.scatter(vcesat["I_C"], vcesat["V_CESAT"])
plt.xscale("log")
plt.xlabel("$I_C$")
plt.ylabel("$V_{CE,SAT}$")
plt.ylim(0, 0.6)
plt.title("$V_{CE,SAT}$ vs $I_C$")
# %%
