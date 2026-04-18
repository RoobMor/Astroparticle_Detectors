# coding: utf-8
# @author: R.M.

from Reduced_Mass import mu
from constants import m_n, N0, rho_chi, v_0, c1, c2, v_light
from numpy import exp, sqrt, pi, array, logspace
from Form_Factor import form_factor
import matplotlib.pyplot as plt


# inizio considerando una sigma_n arbitraria
sigma_n = 1 * 1e-46 # fm^-2

#normalizzo la velocità
v_0 /= v_light

def E_0(m_chi):
    return 1/2 * m_chi * v_0**2

def r(m_chi, A):
    m_t = m_n * A
    return 4 * mu(m_chi, A)**2 / (m_chi * m_t)

def sigma_0(m_chi, A):
    return sigma_n * (mu(m_chi, A) * A / m_n)**2

def R_0(m_chi, A):
    return 2/sqrt(pi) * N0/A * rho_chi/m_chi * sigma_0(m_chi, A) * v_0

def diff_rate(m_chi, A, E_recoil):
    E0r = E_0(m_chi) * r(m_chi, A)
    return R_0(m_chi, A) * form_factor(A, E_recoil)**2 * c1/(E0r) * exp(-c2 * E_recoil / E0r)


# =============================================================================
#  Plot 1
# =============================================================================

plt.xscale("log")
plt.yscale("log")

plt.xlabel(r"$E_{R}$ [keV]")
plt.ylabel(r"$dR/dE_R$")

plt.grid(True)

A      = 40
E_r    = logspace(-1, 2, 200)
Masse_chi  = 1e6 * array([1, 10, 100])
Colors = array(["#ffe81c", "#ff8f05", "#ff3705"])

for m_chi, c in zip(Masse_chi, Colors):
    plt.plot(E_r, diff_rate(m_chi, A, E_r),
             color = c,
             label = r"$m_{\chi}$ =" + f"{m_chi/1e6} GeV")

plt.ylim(1e-60,2*diff_rate(Masse_chi[0], A, 1e-1))
    

plt.legend()
plt.show()