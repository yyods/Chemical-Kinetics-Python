import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define rate constants as functions of temperature T
def rate_constants(T):
    k1 = 4.26e16 * np.exp(-44579 / T)
    k2 = 1.65e9 * (T / 298) ** 4.25 * np.exp(-3890 / T)
    k3 = 8.85e12 * np.exp(-19469 / T)
    k4 = 1.71e12 * (T / 298) ** 2.32 * np.exp(-3414 / T)
    k5 = 1.15e13  # constant
    k6 = 1.45e12  # constant
    alpha = k1 * (3 * k5 + 2 * k6) / (k5 + k6)
    beta = k3 * np.sqrt(k1 / (k5 + k6))
    return k1, k2, k3, k4, k5, k6, alpha, beta

# Define a function to calculate [C2H6](t) from the analytical solution
def ethane_concentration(t, C1_0, alpha, beta):
    # The analytical solution obtained by substituting u(t) = sqrt(C(t))
    # into the transformed ODE leads to:
    return ((-beta/alpha) + np.exp(-alpha/2 * t)*np.sqrt(C1_0) + np.exp(-alpha/2 * t)*(beta/alpha))**2

# Function to perform numerical integration (trapezoidal rule)
def integrate_species(t_vals, func, *args):
    return np.array([np.trapz(func(t_vals[:i], *args), t_vals[:i]) for i in range(1, len(t_vals)+1)])

# Simulation parameters
T = 1100  # Temperature in K
k1, k2, k3, k4, k5, k6, alpha, beta = rate_constants(T)
p = 2e5   # Pressure in Pa
R = 8.3154  # Gas constant in J/(mol K)
C1_0 = p / (R * T * 1e6)  # Adjust units if necessary

time_eval = np.linspace(0, 0.5, 200)  # Time in seconds
C1 = ethane_concentration(time_eval, C1_0, alpha, beta)
C3 = 2 * k1 * integrate_species(time_eval, ethane_concentration, C1_0, alpha, beta)
# For ethylene and hydrogen, we assume proper expressions based on Equations (\ref{eq:C2H4_int}) and (\ref{eq:H2_int})
C6 = beta * integrate_species(time_eval,
          lambda u, C1_0, alpha, beta: np.sqrt(ethane_concentration(u, C1_0, alpha, beta))
          + (k6 / (k3**2)) * ethane_concentration(u, C1_0, alpha, beta),
          C1_0, alpha, beta)
C7 = beta * integrate_species(time_eval,
          lambda u, C1_0, alpha, beta: np.sqrt(ethane_concentration(u, C1_0, alpha, beta)),
          C1_0, alpha, beta)
C8 = (beta**2 * k5 / k3**2) * integrate_species(time_eval, ethane_concentration, C1_0, alpha, beta)

# Plot the results
fig, ax1 = plt.subplots(figsize=(8, 6))
ax2 = ax1.twinx()

ax1.set_ylim([0, 2.5e-5])  # Left y-axis: Ethane, Ethylene, Hydrogen
ax2.set_ylim([0, 1e-6])    # Right y-axis: Methane, Butane

ax1.plot(time_eval, C1, label='Ethane(t)', color='black', linestyle='-')
ax1.plot(time_eval, C6, label='Ethylene(t)', color='black', linestyle='dashed')
ax1.plot(time_eval, C7, label='Hydrogen(t)', color='red', linestyle='dashdot')

ax2.plot(time_eval, C3, label='Methane(t)', color='blue', linestyle='dotted')
ax2.plot(time_eval, C8, label='Butane(t)', color='green', linestyle=':')

ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Concentration (Ethane, Ethylene, Hydrogen)')
ax2.set_ylabel('Concentration (Methane, Butane)')
ax1.set_title('Kinetic Curves at T = 1100 K')

ax1.legend(loc='upper right')
ax2.legend(loc='upper left')
ax1.grid()
plt.show()