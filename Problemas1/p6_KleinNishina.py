import numpy as np
import matplotlib.pyplot as plt

# Constantes
m_p = 938.272  # masa del protón en MeV/c^2
energies = [1, 10, 100, 1000]  # energías del fotón incidente en MeV
theta = np.linspace(0, 2*np.pi, 500)  # ángulo en radianes
theta_deg = np.degrees(theta)       # ángulo en grados para la gráfica cartesiana

# Funciones físicas
def E_prime(E, theta):
    return E / (1 + (E / m_p) * (1 - np.cos(theta)))

def klein_nishina(E, theta):
    Ep = E_prime(E, theta)
    ratio = Ep / E
    return ratio**2 * (ratio + 1/ratio - np.sin(theta)**2)

# Crear figura con dos subgráficas
fig = plt.figure(figsize=(16, 6))

# Gráfica 1: Polar
ax1 = fig.add_subplot(1, 2, 1, polar=True)
for E in energies:
    d_sigma = klein_nishina(E, theta)
    ax1.plot(theta, d_sigma, label=f'{E} MeV')

ax1.set_title('Gráfica polar\nKlein-Nishina γ–protón', va='bottom', fontsize=13)
ax1.set_theta_zero_location('N')
ax1.set_theta_direction(-1)
ax1.set_thetagrids(range(0, 360, 30), labels=[f'{d}°' for d in range(0, 360, 30)])
ax1.set_rgrids([0.5, 1.0, 1.5, 2.0], angle=22.5)
ax1.set_ylim(0, 2.0)
ax1.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))

# Gráfica 2: Cartesiana
ax2 = fig.add_subplot(1, 2, 2)
for E in energies:
    d_sigma = klein_nishina(E, theta)
    ax2.plot(theta_deg, d_sigma, label=f'{E} MeV')

ax2.set_title('Gráfica cartesiana\nSección eficaz vs θ', fontsize=13)
ax2.set_xlabel('Ángulo de dispersión θ (grados)', fontsize=11)
ax2.set_ylabel('Sección eficaz diferencial', fontsize=11)
ax2.grid(True)
ax2.legend()

# Guardar figura combinada
plt.tight_layout()
plt.savefig("graf_polar_cart_KN.png")



