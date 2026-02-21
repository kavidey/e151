# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
# %%
av = pd.read_csv("data/scope_5.csv", header=1)

t = av['second']
vin = av['Volt']
vout = av['Volt.1']
# %%
def wave(t, freq, amplitude, phase, offset):
    return np.sin(t * freq*2*np.pi + phase) * amplitude + offset
# %%
fit_in = curve_fit(wave, av['second'], av['Volt'], p0=[50e3, 0.01, 0, 0])
fit_out = curve_fit(wave, av['second'], av['Volt.1'], p0=[50e3, 1.0, 0, 0])

plt.subplot(2,1,1)
plt.plot(av['second'], av['Volt'])
plt.plot(av['second'], wave(av['second'], *fit_in[0]))

plt.subplot(2,1,2)
plt.plot(av['second'], av['Volt.1'])
plt.plot(av['second'], wave(av['second'], *fit_out[0]))

# print(fit_out[0][0], fit_in[0][0])
print(fit_out[0][1], fit_in[0][1])
fit_out[0][1]/fit_in[0][1]
# %%
r_in = pd.read_csv("data/scope_15.csv", header=1)

t = r_in['second']
vin_before = r_in['Volt.1']
vin_after = r_in['Volt']
# %%
fit_before = curve_fit(wave, t, vin_before, p0=[100e3, 0.01, 0, 0])
fit_after = curve_fit(wave, t, vin_after, p0=[100e3, 1.0, 0, 0])

plt.subplot(2,1,1)
plt.plot(t, vin_before)
plt.plot(t, wave(t, *fit_before[0]))

plt.subplot(2,1,2)
plt.plot(t, vin_after)
plt.plot(t, wave(t, *fit_after[0]))

# print(fit_out[0][0], fit_in[0][0])
print(fit_before[0][1] * 1000, fit_after[0][1]*1000)
# %%
r_out = pd.read_csv("data/scope_9.csv", header=1)

t = r_out['second']
vin = r_out['Volt']
vout = r_out['Volt.1']
# %%
fit_in = curve_fit(wave, t, vin, p0=[50e3, 0.01, 0, 0])
fit_out = curve_fit(wave, t, vout, p0=[50e3, 1.0, 0, 0])

plt.subplot(2,1,1)
plt.plot(t, vin)
plt.plot(t, wave(t, *fit_in[0]))

plt.subplot(2,1,2)
plt.plot(t, vout)
plt.plot(t, wave(t, *fit_out[0]))

# print(fit_out[0][0], fit_in[0][0])
print(fit_in[0][1] * 1000, fit_out[0][1]*1000)
# %%
