import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define rate constants
k1 = 1.0    # Forward rate constant for A -> B
k2 = 0.8    # Reverse rate constant for B -> A
k3 = 0.01   # Slow rate constant for B -> C

# Initial concentrations (in M)
A0 = 0.1
B0 = 0.0
C0 = 0.0

# Define the system of ODEs
def reaction_system(t, y):
    A, B, C = y
    dA_dt = -k1 * A + k2 * B
    dB_dt = k1 * A - (k2 + k3) * B
    dC_dt = k3 * B
    return [dA_dt, dB_dt, dC_dt]

# Time span for simulation
t_span = (0, 15)
t_eval = np.linspace(0, 15, 200)

# Solve the ODEs
solution = solve_ivp(reaction_system, t_span, [A0, B0, C0], t_eval=t_eval)
t = solution.t
A_sol, B_sol, C_sol = solution.y

# Compute the quasi-equilibrium approximation for [B]
B_quasi = (k1 / k2) * A_sol

# Plot the results
plt.figure(figsize=(8, 5))
plt.plot(t, B_sol, label=r'Numerical $B(t)$', linestyle='-', linewidth=2)
plt.plot(t, B_quasi, label=r'Quasi-Equilibrium $\frac{k_1}{k_2}A(t)$', linestyle='dotted', linewidth=2)
plt.xlabel('Time (s)')
plt.ylabel('Concentration (M)')
plt.title('Verification of the Quasi-Equilibrium Approximation')
plt.legend()
plt.grid(True)
plt.show()