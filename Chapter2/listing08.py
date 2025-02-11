from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# Define rate equations for parallel reactions
def parallel_reactions(y, t, k1, k2):
    A, B, C = y
    dA_dt = -(k1 + k2) * A
    dB_dt = k1 * A
    dC_dt = k2 * A
    return [dA_dt, dB_dt, dC_dt]

# Initial conditions
A0, B0, C0 = 1.0, 0.0, 0.0  # Initially, only A is present
k1, k2 = 0.81, 0.27         # Rate constants for the competing pathways
time = np.linspace(0, 10, 100)  # Time range for the simulation

# Solve the system numerically
y0 = [A0, B0, C0]
solution = odeint(parallel_reactions, y0, time, args=(k1, k2))

# Extract the concentration profiles
A_conc = solution[:, 0]
B_conc = solution[:, 1]
C_conc = solution[:, 2]

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(time, A_conc, label=r'$C_A(t)$', linestyle='solid')
plt.plot(time, B_conc, label=r'$C_B(t)$', linestyle='dotted')
plt.plot(time, C_conc, label=r'$C_C(t)$', linestyle='dashed')
plt.xlabel('Time (t)')
plt.ylabel('Concentration')
plt.title('Parallel Reaction Kinetics')
plt.legend()
plt.grid(True)
plt.show()