# %%
%matplotlib widget

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
# %%
a_10x = pd.read_csv("data/scope_4.csv", index_col=None)
a_10x = a_10x.apply(pd.to_numeric, errors='coerce')
a_10x = a_10x.dropna()
a_10x_time = a_10x[np.abs(a_10x['x-axis']) < 1]
a_10x_freq = a_10x[np.abs(a_10x['x-axis']) > 1]

b_10x = pd.read_csv("data/scope_5.csv", index_col=None)
b_10x = b_10x.apply(pd.to_numeric, errors='coerce')
b_10x = b_10x.dropna()
b_10x_time = b_10x[np.abs(b_10x['x-axis']) < 1]
b_10x_freq = b_10x[np.abs(b_10x['x-axis']) > 1]

a_bnc = pd.read_csv("data/scope_6.csv", index_col=None)
a_bnc = a_bnc.apply(pd.to_numeric, errors='coerce')
a_bnc = a_bnc.dropna()
a_bnc_time = a_bnc[np.abs(a_bnc['x-axis']) < 1]
a_bnc_freq = a_bnc[np.abs(a_bnc['x-axis']) > 1]

b_bnc = pd.read_csv("data/scope_7.csv", index_col=None)
b_bnc = b_bnc.apply(pd.to_numeric, errors='coerce')
b_bnc = b_bnc.dropna()
b_bnc_time = b_bnc[np.abs(b_bnc['x-axis']) < 1]
b_bnc_freq = b_bnc[np.abs(b_bnc['x-axis']) > 1]
# %%
a_10x_peaks, _ = find_peaks(a_10x_freq['1'], height=-40, distance=500)
a_10x_peaks = pd.DataFrame({"freq": a_10x_freq['x-axis'].iloc[a_10x_peaks], "height": a_10x_freq['1'].iloc[a_10x_peaks]})
print(a_10x_peaks)

b_10x_peaks, _ = find_peaks(b_10x_freq['1'], height=-40, distance=500)
b_10x_peaks = pd.DataFrame({"freq": b_10x_freq['x-axis'].iloc[b_10x_peaks], "height": b_10x_freq['1'].iloc[b_10x_peaks]})
print(b_10x_peaks)

a_bnc_peaks, _ = find_peaks(a_bnc_freq['1'], height=-40, distance=500)
a_bnc_peaks = pd.DataFrame({"freq": a_bnc_freq['x-axis'].iloc[a_bnc_peaks], "height": a_bnc_freq['1'].iloc[a_bnc_peaks]})
print(a_bnc_peaks)

b_bnc_peaks, _ = find_peaks(b_bnc_freq['1'], height=-40, distance=500)
b_bnc_peaks = pd.DataFrame({"freq": b_bnc_freq['x-axis'].iloc[b_bnc_peaks], "height": b_bnc_freq['1'].iloc[b_bnc_peaks]})
print(b_bnc_peaks)
# %%
fft = np.fft.fft(a_10x_time['1'].to_numpy(), norm="forward")
fftfreq = np.fft.fftfreq(a_10x_time['1'].to_numpy().shape[0], np.gradient(a_10x_time['x-axis'].to_numpy()).mean())
plt.plot(fftfreq, np.abs(fft))
plt.xlim(0, 110e3)
# %%
fig, axs = plt.subplots(1, 2, sharex=True, figsize=(10,5))
axs[0].set_title("Port A")
axs[0].plot(a_10x_time['x-axis'], a_10x_time['1'], label="10x")
axs[0].plot(a_bnc_time['x-axis'], a_bnc_time['1'], label="BNC")
axs[0].set_xlim(0, 0.0001)
axs[0].legend()

axs[1].set_title("Port B")
axs[1].plot(b_10x_time['x-axis'], b_10x_time['1'])
axs[1].plot(b_bnc_time['x-axis'], b_bnc_time['1'])

plt.tight_layout()
# %%
a_10x_peaks['vpp'] = 10**(a_10x_peaks['height']/20) * 2 * np.sqrt(2)
b_10x_peaks['vpp'] = 10**(b_10x_peaks['height']/20) * 2 * np.sqrt(2)
a_bnc_peaks['vpp'] = 10**(a_bnc_peaks['height']/20) * 2 * np.sqrt(2)
b_bnc_peaks['vpp'] = 10**(b_bnc_peaks['height']/20) * 2 * np.sqrt(2)
# %%
C1 = 15e-12
C2 = 60e-12
# %%
porta = pd.DataFrame(a_10x_peaks['freq'])
porta['G1'] = a_10x_peaks['vpp'].to_numpy()
porta['G2'] = a_bnc_peaks['vpp'].to_numpy()
porta['G'] = porta['G1'] / porta['G2']
porta['Z'] = (1-porta['G']) / (porta['G'] * C1 * porta['freq'] * 2 * np.pi - C2 * porta['freq'] * 2 * np.pi)

plt.scatter(porta['freq'], porta['Z'])
plt.xscale('log')
# %%
portb = pd.DataFrame(a_10x_peaks['freq'])
portb['G1'] = b_10x_peaks['vpp'].to_numpy()
portb['G2'] = b_bnc_peaks['vpp'].to_numpy()
portb['G'] = portb['G1'] / portb['G2']
portb['Z'] = (1-portb['G']) / (portb['G'] * C1 * portb['freq'] * 2 * np.pi - C2 * portb['freq'] * 2 * np.pi)

plt.scatter(portb['freq'], portb['Z'])
plt.xscale('log')
# %%
a_10x_zoomed_in = pd.read_csv("data/scope_16.csv", index_col=None, header=1)
a_bnc_zoomed_in = pd.read_csv("data/scope_17.csv", index_col=None, header=1)
# %%
# plt.title("Port A")
plt.plot(a_10x_zoomed_in['second'], a_10x_zoomed_in['Volt'], label="10x")
plt.plot(a_bnc_zoomed_in['second'], a_bnc_zoomed_in['Volt'], label="BNC")
# plt.xlim(0, 0.0001)
plt.legend()
 # %%
