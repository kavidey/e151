# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

%config InlineBackend.figure_format = 'retina'
# %%
def parallel(Z):
    return 1/np.sum(1/np.array(Z))

def left_right(RL, RR, gm):
    return RL + RR + gm * RL * RR

def pole_zero(s, pz):
    return (s/pz) / (s/pz + 1)

def pole(s, p):
    return (1/(s/p + 1))
# %%
f = np.logspace(2.9,7.5, 1000)
s = 1j * f * 2 * np.pi
# %%
Rc = 3.233e3
Re = 5079
Rs = 50
R1 = 20e3
R2 = 10.130e3
R3 = 11.750e3
R4 = 10.090e3
Ratten = 50

ro = 150e3
beta = 200
Ic = 0.97e-3

phi_th = 26e-3
rpi = beta * phi_th / Ic
gm = Ic / phi_th

Ce = 1e-6
Cb = 1e-6
Cmu = 4e-12
Cpi = 15.29e-12
Cbread = 2.15e-12
Cbnc = 49e-12
Cscope = 10e-12
# %%
OCTC_pi = parallel([Rs, Ratten, rpi, R1, R2]) * (2*Cbread + Cpi + Cbnc)
OCTC_mu = left_right(parallel([Rs, Ratten, R1, R2, rpi]), parallel([rpi, (ro + Rc)/(gm*ro)]), gm) * Cmu
OCTC_3 = parallel([rpi, (ro + Rc)/(gm*ro+1), ro]) * (2*Cbread)
OCTC_mu2 = parallel([Rc, beta * ro]) * Cmu
OCTC_pi2 = parallel([rpi, ro, (ro + Rc)/(gm*ro+1)]) * Cpi
OCTC_scope = Rc * (Cbread + Cscope)

print(OCTC_pi, OCTC_mu, OCTC_3, OCTC_mu2, OCTC_pi2, OCTC_scope)

OCTC = OCTC_pi + OCTC_mu + OCTC_3 + OCTC_mu2 + OCTC_pi2 + OCTC_scope
# %%
SCTC_cb1 = Cb * (parallel([Rs, Ratten]) + parallel([rpi, R1, R2]))
SCTC_cb2 = (Cb+2*Cbread) * (beta * ro)

gm1eff = rpi/(rpi + parallel([Rs, Ratten]) ) * gm
SCTC_ce = Ce * parallel([(ro + parallel([1/gm, rpi])) / (gm1eff * ro + 1), parallel([Rs, Ratten, R1, R2]) + rpi, Re])

SCTC = 1/SCTC_cb1 + 1/SCTC_cb2 + 1/SCTC_ce
# %%
octc_approx = (-gm*Rc) * pole(s, 1/OCTC) * pole_zero(s, SCTC)
octc_db = 20 * np.log10(np.abs(octc_approx))
octc_angle = np.angle(octc_approx)
# %%
cascode = pd.read_csv("data/cascode.csv", sep='	')

fig, axs = plt.subplots(2,1, sharex=True)
axs[0].plot(f, octc_db, label="OCTC")
axs[0].set_xscale("log")
axs[0].scatter(cascode["Frequency [kHz]"]*1e3, cascode["dB"])
axs[0].legend()

axs[1].plot(f, np.unwrap(octc_angle/np.pi*180)+180)
axs[1].scatter(cascode["Frequency [kHz]"]*1e3, cascode["Phase"])

axs[0].set_title("Cascode")
# %%
amp_b = pd.read_csv("data/amp_b.csv", sep='	')

fig, axs = plt.subplots(2,1, sharex=True)
axs[0].set_xscale("log")
axs[0].scatter(amp_b["Frequency [kHz]"]*1e3, amp_b["dB"], label='Amp B')
axs[0].scatter(cascode["Frequency [kHz]"]*1e3, cascode["dB"], label='Cascode')
axs[0].legend()

axs[1].scatter(amp_b["Frequency [kHz]"]*1e3, amp_b["Phase (out to in)"])
axs[1].scatter(cascode["Frequency [kHz]"]*1e3, cascode["Phase"])

axs[0].set_title("Cascode vs Amp B")
# %%
