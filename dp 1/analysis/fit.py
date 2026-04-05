# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
# %%
av = pd.read_csv("data/scope_14.csv", header=1)

t = av['second'][:-10]
vin = av['Volt'][:-10]
vout = av['Volt.1'][:-10]
vout2 = av['Volt.1'][:-10]
# %%
def wave(t, freq, amplitude, phase, offset):
    return np.sin(t * freq*2*np.pi + phase) * amplitude + offset
# %%
fit_in = curve_fit(wave, t, vin, p0=[15e3, 0.01, 0, 0])
fit_out = curve_fit(wave, t, vout, p0=[15e3, 1.0, 0, 0])
fit_out2 = curve_fit(wave, t, vout2, p0=[15e3, 1.0, 0, 0])

plt.subplot(3,1,1)
plt.plot(t, vin)
# plt.plot(t, wave(t, *fit_in[0]))

plt.subplot(3,1,2)
plt.plot(t, vout)
plt.plot(t, wave(t, *fit_out[0]))

plt.subplot(3,1,3)
plt.plot(t, vout2)
plt.plot(t, wave(t, *fit_out2[0]))

# print(fit_out[0][0], fit_in[0][0])
print(fit_out[0][1]*2 * 1000, fit_in[0][1]*2 * 1000, fit_out2[0][1]*2 * 1000)
fit_out[0][1]/fit_in[0][1]
# %%
