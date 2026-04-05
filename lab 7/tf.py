# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

%config InlineBackend.figure_format = 'retina'
# %%
def parallel(Z1,Z2):
    return Z1 * Z2 / (Z1 + Z2)
    # return 1/(1/Z1 + 1/Z2)

def tf(s, gm, Rc, ro, Rs, rpi, Cmu, Cpi):
    Rcllro = parallel(Rc,ro)
    Rsllrpi = parallel(Rs,rpi)
    feed_forward_zero = (1-(s*(Cmu/gm)))
    first_order = Cmu*(Rcllro + Rsllrpi + gm*Rcllro*Rsllrpi) + Cpi*Rsllrpi
    second_order = Cmu*Cpi*Rcllro*Rsllrpi
    print(first_order, second_order)
    return (-gm*Rc) *  feed_forward_zero / (1 + s*first_order + (s**2)*second_order)

def miller(gm, s, Rs, rpi, Rb, Rc, Cmu, Cpi):
    p_mu = 1/(parallel(parallel(Rs, rpi), Rb) * Cmu * (1 + gm*Rc))
    p_pi = 1/(parallel(parallel(Rs, rpi), Rb) * Cpi)
    return (-gm*Rc) * 1/(s/p_mu + 1) * 1/(s/p_pi + 1)
# %%
f = np.logspace(2,7.5, 1000)
s = 1j * f * 2 * np.pi
# %%
gm = 39.97e-3
ro = 150e3
Rc = 3.233e3
Rb = 990
Re = 5079
beta = 200
Ic = 1.039e-3
phi_th = 26e-3
rpi = beta * phi_th / Ic
Rs = 50
rin = 881.6
Cmu = 4e-12
Cpi = 15.29e-12
Cbread = 2.15e-12
# %%
pz1 = 1/(1e-6 * (rin+50)) # 1 uF coupling capacitor sees rin of CE + 50Ω src impedance of fgen
tf_gain = tf(s, gm, Rc, ro, Rs, rpi, Cmu + Cbread, Cpi + Cbread) * (s/pz1)/(s/pz1 + 1)
tf_db = 20 * np.log10(np.abs(tf_gain))
tf_angle = np.angle(tf_gain)

p_bread = 1/(Cbread * Rc)
miller_gain = miller(gm, s, Rs, rpi, Rb, Rc, Cmu + Cbread, Cpi + Cbread) * (s/pz1)/(s/pz1 + 1)

miller_db = 20 * np.log10(np.abs(miller_gain))
miller_angle = np.angle(miller_gain)

amp_a = pd.read_csv("data/amp_a.csv", sep='	')

fig, axs = plt.subplots(2,1, sharex=True)
axs[0].plot(f, tf_db, label="Full TF")
axs[0].plot(f, miller_db, label="Miller")
axs[0].set_xscale("log")
axs[0].scatter(amp_a["Frequency [kHz]"]*1e3, amp_a["dB"])
axs[0].legend()

axs[1].plot(f, np.unwrap(tf_angle/np.pi*180) + 180)
axs[1].plot(f, np.unwrap(miller_angle/np.pi*180) + 180)
axs[1].scatter(amp_a["Frequency [kHz]"]*1e3, amp_a["Phase (out to in)"])

axs[0].set_title("Amp A")
# %%
pz1 = 1/(1e-6 * (rin+50)) # 1 uF coupling capacitor sees rin of CE + 50Ω src impedance of fgen
tf_gain = tf(s, gm, Rc, ro, Rs, rpi, Cmu, Cpi + 2*Cbread) * (s/pz1)/(s/pz1 + 1)
tf_db = 20 * np.log10(np.abs(tf_gain))
tf_angle = np.angle(tf_gain)

p_bread = 1/(Cbread * Rc)
miller_gain = miller(gm, s, Rs, rpi, Rb, Rc, Cmu, 2*Cbread) * (s/pz1)/(s/pz1 + 1) * 1/(s/p_bread + 1)
miller_db = 20 * np.log10(np.abs(miller_gain))
miller_angle = np.angle(miller_gain)

amp_a = pd.read_csv("data/amp_b.csv", sep='	')

fig, axs = plt.subplots(2,1, sharex=True)
axs[0].plot(f, tf_db, label="Full TF")
axs[0].plot(f, miller_db, label="Miller")
axs[0].set_xscale("log")
axs[0].scatter(amp_a["Frequency [kHz]"]*1e3, amp_a["dB"])
axs[0].legend()

axs[1].plot(f, np.unwrap(tf_angle/np.pi*180) + 180)
axs[1].plot(f, np.unwrap(miller_angle/np.pi*180) + 180)
axs[1].scatter(amp_a["Frequency [kHz]"]*1e3, amp_a["Phase (out to in)"])

axs[0].set_title("Amp B")
# %%
