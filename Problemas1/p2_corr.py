import numpy as np

# Datos físicos (CGS)
rho_Ag = 10.49  # g/cm³
molar_m_Ag = 107.87  # g/mol
n_disp_cent = rho_Ag * 6.022e23 / molar_m_Ag  # centros dispersores / cm³

# Constantes y condiciones del problema
Z1, Z2 = 2, 47
T = 5e6 * 1.602e-12  # energía en ergios
s = 0.2  # cm (ancho del detector)
r = 20   # cm (radio del detector)
area = 1  # cm²
thick = 1e-4  # cm (espesor del blanco)
flux = 1e10 * 60  # partículas / (cm²·min)
e2 = 1.44e-13  # erg·cm

def sigma_total(theta_deg, Z1, Z2, T, s, r):
    theta = np.deg2rad(theta_deg)
    theta1 = theta / 2
    theta2 = theta1 + s / (2 * r)  # ya en radianes
    const = 4 * np.pi * ((Z1 * Z2 * e2) / (4 * T))**2  # en cm²
    return const * (1 / np.sin(theta1)**2 - 1 / np.sin(theta2)**2)

def dNdt(sigma, n, A, dx, flux):
    return sigma * n * A * dx * flux

# Cálculo para dos ángulos
theta_1 = 45
theta_2 = 135

sigma1 = sigma_total(theta_1, Z1, Z2, T, s, r)
sigma2 = sigma_total(theta_2, Z1, Z2, T, s, r)

N1 = dNdt(sigma1, n_disp_cent, area, thick, flux)
N2 = dNdt(sigma2, n_disp_cent, area, thick, flux)

print(f"sigma({theta_1}°) = {sigma1:.3e} cm²")
print(f"sigma({theta_2}°) = {sigma2:.3e} cm²")
print(f"dN/dt({theta_1}°) = {N1:.3e} partículas/min")
print(f"dN/dt({theta_2}°) = {N2:.3e} partículas/min")