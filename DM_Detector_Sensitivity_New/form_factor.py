import numpy as np
import math
import matplotlib.pyplot as plt
from constants import *

# UNITA' DI RIFERIMENTO PER I SEGUENTI CONTI: keV, fm

def form_factor(q, r_0, s):
    t = q*r_0
    num = 3 * (np.sin(t) - t*np.cos(t))
    den = t ** 3
    expon = np.exp(-(q*s)**2/2)
    return expon * num / den

def f_squared_calculation(Er, A):
    s      = 0.9                                # fm
    r_n = 1.2 * A**(1/3)                        # raggio nucleare, fm
    r_0 = math.sqrt(5*(r_n**2 / 3 - s**2))      # fm
    r_0 = r_0 * fm_to_1_over_keV                # keV^-1
    m_N = A * m_protone                         # keV
    q   = np.sqrt(2 * Er * m_N)                 # keV

    F = form_factor(q, r_0, s* fm_to_1_over_keV)
    F_squared = F**2

    return F_squared

def f_plotting(E_th, E_max, plotting_flag, A_list = [40, 72, 131]):

    colors = plt.cm.plasma(np.linspace(0.6, 1, len(A_list)))

    fig, ax = plt.subplots()
    ax.set_title(fr"Squared Form Factor")
    ax.set_xlabel(fr"$E_R$ [keV]")
    ax.set_ylabel(fr"$|F|^2$")

    for A, color in zip(A_list, colors):
        Er  = np.logspace(np.log10(E_th), np.log10(E_max), 2000)
        F_squared = f_squared_calculation(Er, A)

        ax.plot(Er, F_squared, label = fr"$A = {A}$", color = color)

    if plotting_flag:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid(True, alpha = 0.3)
        ax.legend()
    else:
        plt.close()

    return F_squared

# f_plotting(1e-1, 1e4, 1)
# plt.show()