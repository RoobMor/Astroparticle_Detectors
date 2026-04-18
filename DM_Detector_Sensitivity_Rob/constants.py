# coding: utf-8
# @author: R.M.

from numpy import pi

# velocità in km/s
v_c    = 238.
v_odot = v_c + 12.
v_plus = 29.8
v_esc  = 544.
v_0 = 220.
v_light = 2.998e5

# parametri orbitali
omega = 2 * pi / 365 # periodo orbitale terrestre, 1/day
theta = pi/3         # inclinazione orbita terrestre - piano galattico

# fattori adimensionali
N0 = 6.022e29       # N. Avogadro
c1 = 0.751
c2 = 0.561

# fattori di conversione
s_to_yr  = 1/(365*24*60*60)
km_to_cm = 1e5
days_to_years = 1/365
fm_to_1_over_keV = 1/197.3269631 * 1e-3

# massa nucleoni, keV
m_n = 939 * 1e3

# skin factor, fm
s = 0.9 * fm_to_1_over_keV

# densità di massa della particella chi, keV/fm^3
# rho_chi = 0.3 GeV/cm^3
rho_chi = 0.3 * 1e6 * 1e-39
