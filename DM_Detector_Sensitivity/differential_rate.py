import numpy as np
import math
import matplotlib.pyplot as plt
from constants import *
from form_factor import f_squared_calculation, f_plotting
from scipy.integrate import trapezoid


def dR_dE_calculation(Er, A, F2, m_chi, mu, E0, r):
    # 1 / (ton yr GeV)
    return 2/math.sqrt(pi) * NA/A * rho_chi/m_chi * sigma_0 * (mu/(m_protone * 1e-6) * A)**2 * v0 * F2 * c1/(E0 * r) * np.exp(-c2*Er/(E0*r)) * 1e38 / s_to_yr


def sigma_exclusion_calculation(m_chi, A, E_th, MT):
    m_N   = A * m_protone * 1e-6                                        # GeV
    mu    = m_chi * m_N / (m_chi + m_N)                                 # GeV
    r     = 4 * mu**2 / (m_chi * m_N)                                   # -

    E0    = 0.5 * m_chi * v0**2                                         # GeV
    E_max = 2 * mu**2 * v_max**2 / m_N                                  # GeV
    E_min = E_th * 1e-6                                                 # GeV
    Er    = np.logspace(np.log10(E_min), np.log10(E_max[0]), 1000)      # GeV

    # Nel fattore di forma, Er deve entrare in keV
    F2    = f_squared_calculation(Er*1e6, A)

    n_expected = MT * trapezoid(efficiency * dR_dE_calculation(Er, A, F2, m_chi, mu, E0, r), Er, axis=0) # events
    # ATTENZIONE!
    # trapezoid() restituisce un singolo valore -> (valore)
    # quad() restituisce una tupla -> (valore, errore_stimato)

    return 2.3 * sigma_0 / n_expected # cm^2


# UNITA' DI RIFERIMENTO PER I SEGUENTI CONTI: keV, fm

# Listed Parameters
A_list    = [40, 72, 131]          # Numero di massa
MT_list   = [2, 20, 200]           # ton yr
E_th_list = [2, 10, 50]            # keV


# Fixed Parameters
A_    = 40              # Argon
M_    = 100             # ton
T_    = 1               # yr
E_th_ = 1               # keV

# Constants
sigma_0    = 1                                       # cm^2
efficiency = 1                                       # -
m_chi = np.logspace(0, 4, 1000)[None, :]             # GeV


colors = plt.cm.plasma(np.linspace(0.6, 1, len(A_list)))


f_plotting(E_th_, 1e4, 1, A_list = A_list)


fig1, ax1 = plt.subplots()

for A, color in zip(A_list, colors):
    sigma_exclusion = sigma_exclusion_calculation(m_chi, A, E_th_, M_ * T_) # cm^2
    ax1.plot(m_chi[0], sigma_exclusion, label = fr"$A = {A}$", color = color)

ax1.set_title(fr"Exclusion Limit: Vanilla Plot. $\ MT = {M_ * T_}$ ton yr, $\ E_{{th}} = {E_th_}$ keV.")
ax1.set_xlabel(fr"$m_\chi$ [GeV]")
ax1.set_ylabel(fr"$\sigma$ [cm$^2$]")
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.grid(True, alpha = 0.3)
ax1.legend()


fig2, ax2 = plt.subplots()

for MT, color in zip(MT_list, colors):
    sigma_exclusion = sigma_exclusion_calculation(m_chi, A_, E_th_, MT) # cm^2
    ax2.plot(m_chi[0], sigma_exclusion, label = fr"$MT = {MT}$ ton yr", color = color)

ax2.set_title(fr"Exclusion Limit: Vanilla Plot. $\ A = {A_}$ (Argon), $\ E_{{th}} = {E_th_}$ keV.")
ax2.set_xlabel(fr"$m_\chi$ [GeV]")
ax2.set_ylabel(fr"$\sigma$ [cm$^2$]")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.grid(True, alpha = 0.3)
ax2.legend()


fig3, ax3 = plt.subplots()

for E_th, color in zip(E_th_list, colors):
    sigma_exclusion = sigma_exclusion_calculation(m_chi, A_, E_th, M_ * T_) # cm^2
    ax3.plot(m_chi[0], sigma_exclusion, label = fr"$E_{{th}} = {E_th}$ keV", color = color)

ax3.set_title(fr"Exclusion Limit: Vanilla Plot. $\ MT = {M_ * T_}$ ton yr, $\ A = {A_}$ (Argon).")
ax3.set_xlabel(fr"$m_\chi$ [GeV]")
ax3.set_ylabel(fr"$\sigma$ [cm$^2$]")
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.grid(True, alpha = 0.3)
ax3.legend()

plt.show()