from constants import s, m_n, fm_to_1_over_keV
from numpy import sin, cos, exp, power, sqrt


def J1(z):
    return (sin(z) - z * cos(z)) / z**2

def form_factor(A, E_r):
    r_n = 1.2 * power(A, 1/3)
    r_0 = sqrt(5*(r_n**2 / 3 - s**2))
    m_N = m_n * A
    q = sqrt(2 * E_r * m_N) * fm_to_1_over_keV
    
    return 3 * J1(q*r_0) / (q*r_0) * exp(-(q*s)**2 / 2)

def mu(m_chi, A):
    m_t = A * m_n
    return (m_chi * m_t) / (m_chi + m_t)