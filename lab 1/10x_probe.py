# %%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.signal import find_peaks, hilbert, sosfilt, butter

%config InlineBackend.figure_format = 'retina'
# %%
def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
# %%
bnc_in = pd.read_csv("data/scope_8.csv", header=1)
bnc_out = pd.read_csv("data/scope_9.csv", header=1)

in_t = bnc_in['second'].to_numpy()
in_v = bnc_in['Volt'].to_numpy()
in_v = in_v[(in_t>0) & (in_t<1e-3)]
in_t = in_t[(in_t>0) & (in_t<1e-3)]


out_t = bnc_out['second'].to_numpy()
out_v = bnc_out['Volt'].to_numpy()
out_v = out_v[(out_t>0) & (out_t<1e-3)]
out_t = out_t[(out_t>0) & (out_t<1e-3)]
# %%
plt.plot(in_t+0.000002, in_v)
plt.plot(out_t, out_v)
# %%
# f = np.logspace(np.log10(80e3), np.log10(1.5e6), in_t.shape[0])
f = np.logspace(np.log10(100e3), np.log10(1.5e6), in_t.shape[0])
# f = np.linspace(80e3, 1.5e6, in_t.shape[0])
# %%
in_h = hilbert(in_v)
out_h = hilbert(out_v)

in_amp = np.abs(in_h)
# in_phase = np.unwrap(np.angle(in_h))
in_phase = np.angle(in_h)
out_amp = np.abs(out_h)
# out_phase = np.unwrap(np.angle(out_h))
out_phase = np.angle(out_h)
# %%
plt.plot(f, out_v)
plt.plot(f, out_amp, alpha=0.5)
plt.axhline(0.5 * 0.7, c='orange')
plt.axvline(1.05, c='orange')
# %%
plt.plot(f, in_amp)
plt.plot(f, out_amp)
# %%
sos = butter(10, 2, 'lp', fs=1000, output='sos')
# %%
fig, axs = plt.subplots(2,1, sharex=True)
plt.suptitle("10x Probe Bode Plot")

axs[0].set_xscale('log')
axs[0].set_xlim(2.1e5, 1.5e6)
axs[0].plot(f, sosfilt(sos, 20 * np.log10(out_amp/0.5)), label='raw')
axs[0].plot(f, 20 * np.log10(out_amp/0.5), alpha=0.5, label='filtered')
axs[0].axhline(-3, linestyle='--', c='r', label="-3db")
axs[0].axvline(1.05e6, c='k', label="Measured $f_c$")
axs[0].legend()

phase = np.mod(out_phase-in_phase, 2*np.pi) * 180 / np.pi - 270
# phase = (out_phase-in_phase) * 180 / np.pi + 76
axs[1].plot(f, sosfilt(sos, phase))
axs[1].plot(f, phase, alpha=0.5)
axs[1].set_ylim(-100, 10)
# axs[1].axhline(-45, linestyle='--', c='r')
axs[1].axvline(1.05e6, c='k')
axs[1].xaxis.set_minor_formatter(matplotlib.ticker.LogFormatterSciNotation(minor_thresholds=(2, 0.4)))
# axs[1].xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter(minor_thresholds=(2, 0.4)))

plt.tight_layout()
# %%
print(f[find_nearest_idx(sosfilt(sos, 20 * np.log10(out_amp/0.5)), -3)])
# %%
