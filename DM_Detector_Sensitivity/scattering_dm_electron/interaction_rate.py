import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import erf
from constants import *

#####################################################################
#                                                                   #
# BLOCCO 0 - DICHIARAZIONE COSTANTI                                 #
#                                                                   #
#####################################################################

# Dichiaro di esporre il prefattore in ng, cm, h
# COSTANTI CARATTERISTICHE (ref. QROCODILE)
omega_threshold = 0.11          # eV
rho_WSi         = 9.3 * 1e9     # ng/cm^3
efficiency      = 0.98          # -
M               = 1.67          # ng
T               = 415.15        # hours


# Fermi parameters
E_f = 7                                 # eV
v_f = np.sqrt(2 * E_f / m_elettrone)    # adim
k_f = m_elettrone * v_f                 # eV


# Massa DM di riferimento
m_chi = 1e4                 # eV = 10 keV


# Massa del mediatore (2 regimi estremali)
mediator_mass = "HM" # Heavy Mediator
# mediator_mass = "LM" # Light Mediator


# Momento trasferito
q = 10                      # eV
q0 = m_elettrone/alpha      # eV


# Sezione d'urto arbitraria di riferimento / normalizzazione
sigma_e_bar = 1e-39         # cm^2, vedi 2101.08263, pag.11, didascalia fig.S4


# Range energetico dell'esperimento
omega = np.logspace(-6,2,1000) # eV


# Range di velocità della particella, standard halo model + cinematica
v_E_modulus = np.sqrt(v_E[0]**2 + v_E[1]**2 + v_E[2]**2)        # adim
v_min       = np.sqrt(2*omega_threshold/m_chi)                  # adim
v_max       = v_esc + v_E_modulus                               # adim

print(v_min)
print(v_max)


#####################################################################
#                                                                   #
# BLOCCO 1 - FUNZIONE DINAMICA DI STRUTTURA                         #
#                                                                   #
#####################################################################


def epsilon(q, omega):
    omega_p = (v_f * e) / (np.sqrt(3) * np.pi) * (2 * E_f * m_elettrone**3)**0.25
    Gamma_p = 0.1 * E_f
    Q_plus  = q/(2*k_f) + (omega + 1j*Gamma_p)/(q*v_f)
    Q_minus = q/(2*k_f) - (omega + 1j*Gamma_p)/(q*v_f)
    prefactor = (3*omega_p**2)/(q**2 * v_f**2)
    add1 = k_f/(4*q) * (1 - Q_minus**2) * np.log((Q_minus+1)/(Q_minus-1))
    add2 = k_f/(4*q) * (1 - Q_plus**2)  * np.log((Q_plus +1)/(Q_plus -1))
    return 1 + prefactor * (1/2 + add1 + add2)


def S(q, omega):
    W = (-1/epsilon(q, omega)).imag # Parte immaginaria del reciproco della funzione dielettrica di struttura
    return 2 * (q / e)**2 * W


fig, ax = plt.subplots()
ax.set_title(fr"Structure Function vs Energy")
ax.set_xlabel(fr"$\omega$ [eV]")
# ax.set_ylabel(fr"$S \ [eV^2]$")
ax.set_ylabel(fr"$W$")
ax.plot(omega, (2 * (q / e)**2)**-1 * S(q, omega))
ax.set_xscale("log")
ax.set_yscale("log")
ax.grid(True, alpha = 0.3)
# ax.legend()


#####################################################################
#                                                                   #
# BLOCCO 2 - MAXWELL-BOLTZMANN INTEGRATA                            #
#                                                                   #
#####################################################################


def domain(value):
    if value >= v_min and value <= v_max:
        return 1
    else:
        return 0


def maxwell_boltzmann(v, cosTheta):
    z = v_esc/v0
    N_E = (erf(z) - 2/np.sqrt(pi) * z * np.exp(-z**2)) * (np.sqrt(pi) * v0)**3
    v_modulus2 = v**2 + v_E_modulus**2 + 2 * v * v_E_modulus * cosTheta
    f_MB = 1/N_E * np.exp(-(v_modulus2)/(v0**2)) * ((v_esc - np.sqrt(v_modulus2)) >= 0)
    return f_MB


def f_tilde(v):
    val, err = quad(maxwell_boltzmann, -1, 1, args=(v,))
    return val


#####################################################################
#                                                                   #
# BLOCCO 3 - RATE DI INTERAZIONE                                    #
# dv dcosTheta dq                                                   #
#####################################################################


def form_factor(q):
    if   mediator_mass == "HM": # Heavy Mediator
        return 1
    elif mediator_mass == "LM": # Light Mediator
        return (q0/q)**2


# Primo integrale: dq
def integrale_dq(v, cosTheta):
    q_upper     = m_chi*v*cosTheta + np.sqrt((m_chi*v*cosTheta)**2 - 2*omega_threshold*m_chi)
    q_lower     = m_chi*v*cosTheta - np.sqrt((m_chi*v*cosTheta)**2 - 2*omega_threshold*m_chi)

    # Nel momento in cui integro il rate differenziale, scatta la delta in omega_p.
    # DEVO GESTIRE QUESTO ASPETTO SE VOGLIO IL RATE INTEGRATO
    omega_p     = lambda q: q * v * cosTheta - q**2 / (2*m_chi)

    integranda  = lambda q: form_factor(q)**2 * f_tilde(v) * v**2 * q**2 * S(q, omega_p(q))

    val, err    = quad(integranda, q_lower, q_upper)
    return val


# Secondo integrale: dcosTheta
def integrale_dcos(v):
    cosTheta_min    = v_min/v
    integranda      = lambda cosTheta: integrale_dq(v, cosTheta)
    val, err        = quad(integranda, cosTheta_min, 1)
    return val


# Terzo integrale: dv
def integrale_dv():
    integranda  = lambda v: integrale_dcos(v)
    val, err    = quad(integranda, v_min, v_max)
    return val



n_chi       = rho_chi/m_chi
mu          = m_chi * m_elettrone / (m_chi + m_elettrone)
prefactor   = pi/2 * n_chi * sigma_e_bar / (mu**2 * rho_WSi)

Gamma = prefactor * integrale_dv()
print(Gamma)
# dR_domega = [prefactor * integrale_dv(w) for w in omega]


# fig, ax = plt.subplots()
# ax.set_title(fr"Differential Rate vs Energy")
# ax.set_xlabel(fr"$\omega$ [eV]")
# ax.set_ylabel(fr"$\diff{{R}}{{\omega}}$ [ng hours]")
# ax.plot(omega, dR_domega)
# ax.set_xscale("log")
# ax.set_yscale("log")
# ax.grid(True, alpha = 0.3)

plt.show()