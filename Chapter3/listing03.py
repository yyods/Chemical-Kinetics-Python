import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define rate constants and initial concentration
k1, k2, k3 = 0.15, 0.07, 10
A0 = 0.1

# Define the system of ODEs
def reaction_system(t, y):
    A, B = y
    dA_dt = -k1 * A
    dB_dt = k1 * A - (k2 + k3) * B
    return [dA_dt, dB_dt]

# Solve the ODEs over the interval from 0 to 2 seconds
t_span = [0, 2]
t_eval = np.linspace(0, 2, 100)
sol = solve_ivp(reaction_system, t_span, [A0, 0], t_eval=t_eval)

# Extract the concentration profiles for A and B
A, B = sol.y
t = sol.t

# Compute the steady-state prediction for the rate of product formation
rP_steady = (k1 * k3 / (k2 + k3)) * A

# Plot the numerical and steady-state predicted rates
plt.figure(figsize=(8, 5))
plt.plot(t, k3 * B, label=r'$k_3 B(t)$', color='black')
plt.plot(t, rP_steady, '--', label=r'$\frac{k_1 k_3}{k_2 + k_3} A(t)$', color='black')
plt.xlabel('Time (s)')
plt.ylabel('Rate (mol/L/s)')
plt.legend()
plt.grid()
plt.show()