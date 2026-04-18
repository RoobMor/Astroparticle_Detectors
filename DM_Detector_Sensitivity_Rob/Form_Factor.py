# coding: utf-8
# @author: R.M.
from constants import s, m_n, fm_to_1_over_keV
import numpy as np
from numpy import sin, cos, exp, power, sqrt
import matplotlib.pyplot as plt

def J1(z):
    return (sin(z) - z * cos(z)) / z**2

def form_factor(A, E_r):
    r_n = 1.2 * power(A, 1/3)
    r_0 = sqrt(5*(r_n**2 / 3 - s**2))
    m_N = m_n * A
    q = sqrt(2 * E_r * m_N) * fm_to_1_over_keV
    
    return 3 * J1(q*r_0) / (q*r_0) * exp(-(q*s)**2 / 2)
 
    
# =============================================================================
#  Plot 1
# =============================================================================
plt.xscale("log")

plt.xlabel(r"$E_R$ [keV]")
plt.ylabel(r"$F^2(E_R)$")

plt.grid(True)

A = np.array([40, 72, 131])   # Numeri di Massa per Argon, Germanio, Xenon
Colors = np.array(["#ffe81c", "#ff8f05", "#ff3705"])
Elementi = np.array(["Argon", "Germanio", "Xenon"])

E_r = np.logspace(-1, 4, 200)
for a, c, e in zip(A, Colors, Elementi):
    max_len = max(len(stringa) for stringa in Elementi)
    padding = '~' * (max_len - len(e) + 1)
    plt.plot(E_r, form_factor(a, E_r)**2,
             color = c,
             label = f"$A_{{{e}}}{padding}$ = {a}")
   
plt.legend()
plt.show()


# =============================================================================
# Plot 2
# =============================================================================
plt.xscale("log")
plt.yscale("log")

plt.xlabel(r"$E_R$ [keV]")
plt.ylabel(r"$F^2(E_R)$")

plt.grid(True)

A = np.array([40, 72, 131])   # Numeri di Massa per Argon, Germanio, Xenon
Colors = np.array(["#ffe81c", "#ff8f05", "#ff3705"])
Elementi = np.array(["Argon", "Germanio", "Xenon"])

E_r = np.linspace(1e-1, 1e4, 5000)
for a, c, e in zip(A, Colors, Elementi):
    max_len = max(len(stringa) for stringa in Elementi)
    padding = '~' * (max_len - len(e) + 1)
    plt.plot(E_r, form_factor(a, E_r)**2,
             color = c,
             label = f"$A_{{{e}}}{padding}$ = {a}")

plt.legend()
plt.show()