from numpy import pi

fm_to_1_over_keV = 1/197.3269631 * 1e-3
s_to_yr          = 1/(365*24*60*60)

m_protone = 938.272 * 1e3 # keV

# Numero di Avogadro
# NA = 6.022e26       # kg^-1
NA = 6.022e29         # ton^-1

c1 = 0.751
c2 = 0.561

# Densità di Materia Oscura
# rho_chi = 0.3 GeV/cm^3
rho_chi = 0.3 * 1e6 * 1e-39 # keV / fm^3

# Per ottenere le velocità in unità naturali, devono essere divise per c
c = 299792.458 # km/s
v0 = 220/c # km/s
v_max = (544 + 220 + 12 + 30)/c # km/s