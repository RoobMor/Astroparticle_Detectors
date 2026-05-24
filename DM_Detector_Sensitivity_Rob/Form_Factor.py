# coding: utf-8
# @author: R.M.

from functions import form_factor
import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
#  Plot 1
# =============================================================================
plt.figure(1)
plt.xscale("log")

plt.title("Form Factors in LogSpace")
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
plt.figure(2)
plt.xscale("log")
plt.yscale("log")

plt.title("Form Factors in LogLogSpace")
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