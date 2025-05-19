import numpy as np
import matplotlib.pyplot as plt

# Constantes físicas
m_p = 938.272  # masa del protón en MeV/c^2
alpha = 1/137  # constante de estructura fina (adimensional)
hbar_c = 197.327  # MeV·fm (ℏ·c en unidades naturales)

# Radio clásico del protón (en fm)
r_p = alpha * hbar_c / m_p  # resultado en femtómetros
r_p_m = r_p * 1e-15  # convertido a metros

# Energías del fotón incidente (en MeV)
energies = [1, 10, 100, 1000]

# Ángulo theta (en radianes)
theta = np.linspace(0, 2* np.pi, 500)

# Función para calcular E' (energía del fotón dispersado)
def E_prime(E, theta):
    return E / (1 + (E / m_p) * (1 - np.cos(theta)))

# Función para calcular la sección eficaz diferencial (Klein-Nishina para γ-protón)
def klein_nishina(E, theta):
    Ep = E_prime(E, theta)
    ratio = Ep / E
    # Fórmula completa con el radio clásico del protón
    return 0.5 * r_p_m**2 * ratio**2 * (ratio + 1/ratio - np.sin(theta)**2)

# Crear gráfica polar
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, polar=True)

# Graficar para cada energía
for E in energies:
    d_sigma = klein_nishina(E, theta)
    ax.plot(theta, d_sigma, label=f'{E} MeV')

# Configuración de la gráfica
ax.set_title('Distribución angular (Klein-Nishina γ–protón)', va='bottom')
ax.set_theta_zero_location('N')  # 0° arriba
ax.set_theta_direction(-1)       # Sentido horario
ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))


# Guardar figura combinada
plt.savefig('graf_polar_rp_KN.png')


