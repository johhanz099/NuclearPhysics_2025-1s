import numpy as np

# Parámetros físicos
E_gamma = 50        # MeV
m_p = 938.27        # MeV
r_p_cm = 1.535e-16  # radio clásico del protón en cm

# Resultado de theta_max calculado anteriormente
tan_phi_target = np.tan(np.radians(5))  # tangente del ángulo de dispersión del protón

# Función f(theta)
def f(theta):
    return (m_p * np.sin(theta)) / ((E_gamma + m_p) * (1 - np.cos(theta))) - tan_phi_target

# Derivada de f(theta)
def df(theta):
    num = m_p * (np.cos(theta)*(1 - np.cos(theta)) - np.sin(theta)**2)
    den = (E_gamma + m_p) * (1 - np.cos(theta))**2
    return num / den

# Método de Newton-Raphson
def newton_raphson(theta0, tol=1e-6, max_iter=100):
    theta = theta0
    for i in range(max_iter):
        f_val = f(theta)
        df_val = df(theta)
        if df_val == 0:
            raise ValueError("Derivada cero. Newton-Raphson falla.")
        theta_next = theta - f_val / df_val
        if abs(theta_next - theta) < tol:
            return theta_next
        theta = theta_next
    raise ValueError("No convergió en el número máximo de iteraciones.")

# Buscar theta_max
theta0 = np.radians(5)
theta_max = newton_raphson(theta0)

# ---- Cálculo de sigma con la fórmula integrada ----
k = E_gamma / m_p
rp2 = r_p_cm**2
theta_pi = np.pi

def f_KN(theta):
    return 1 + k * (1 - np.cos(theta))

f_theta_max = f_KN(theta_max)
f_pi = f_KN(theta_pi)

# Recalculamos sigma usando la fórmula corregida

# Nuevos términos según la corrección
term1 = 0.5 * (1 / f_theta_max**2 - 1 / f_pi**2)
term2 = ((2 * k + 1) / k**2) * (1 / f_theta_max - 1 / f_pi)
term3 = (1 - (2 * (k + 1)) / k**2) * np.log(f_pi / f_theta_max)
term4 = (1 / k**2) * (f_pi - f_theta_max)

# Sección eficaz corregida
sigma_cm2 = (np.pi * rp2 / k) * (term1 + term2 + term3 + term4)
sigma_cm2

# Resultados
print(np.degrees(theta_max), sigma_cm2)

# n_p 
n_parafina = 0.9 * 6.022e23 * 52 /  352.66
print(n_parafina)

# N_p
N_p = 1e11 * n_parafina * 0.1 * sigma_cm2
print(N_p)
