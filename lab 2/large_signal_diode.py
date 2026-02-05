# %%
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
import scipy.stats
import numpy as np
# %%
df = pd.read_csv("data/large_signal_diode.csv", sep="\t")
phi_th = 0.026 # 26 mV

df
# %%
def diode_eqn(v_D, I_s, n):
    return I_s * (np.exp(v_D/(n*phi_th)) - 1)
# %%
popt, pcov = scipy.optimize.curve_fit(
    diode_eqn,
    df["Vd"],
    df["Id"],
    absolute_sigma=False,
    sigma=10,
    p0=[4.5e-10, 1.56], # starting parameters
    bounds=((0, 1), (1, 10)), # parameter bounds
)
# %%
plt.scatter(df["Vd"], df["Id"])
plt.yscale("log")

vd = np.linspace(np.min(df["Vd"]), np.max(df["Vd"]), 100)
plt.plot(vd, diode_eqn(vd, *popt))
print(f"Is = {popt[0]:.4}")
print(f"n = {popt[1]:.4}")
plt.plot(vd, diode_eqn(vd, 3e-10, 1.6))
# %%
res = scipy.stats.linregress(df["Vd"], np.log(df["Id"]))
# %%
plt.yscale("log")

plt.scatter(df["Vd"], df["Id"])

plt.plot(vd, np.exp(res.intercept + res.slope*vd), c='tab:orange', label='fit')

plt.grid()
plt.xlabel("$v_D$")
plt.ylabel("$i_D$")

plt.legend()

Is = np.exp(res.intercept)
n = 1/(res.slope*phi_th)
print(f"Is = {Is:.4}")
print(f"n = {n:.4}")

# plt.plot(vd, diode_eqn(vd, Is, n))
# %%
