import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define reaction rate constants
k1 = 0.5   # Rate constant for A -> B
k2 = 0.1   # Rate constant for B + C -> P

# Initial concentrations (Reordered for better logical flow)
A0 = 1.0   # Initial concentration of A
B0 = 0.0   # Initially, no B is present
C0 = 0.9   # Initial concentration of C
P0 = 0.0   # Initially, no P is present

# Time range for the simulation
t = np.linspace(0, 20, 200)  # Time from 0 to 20 units

# Define the system of differential equations
def reaction_rates(y, t, k1, k2):
    A, B, C, P = y  # Now the variables follow the same order as the initial conditions
    dA_dt = -k1 * A                  # A -> B
    dB_dt = k1 * A - k2 * B * C       # B is formed and consumed in B + C -> P
    dC_dt = -k2 * B * C               # C is consumed in the reaction
    dP_dt = k2 * B * C                # P is formed from B + C
    return [dA_dt, dB_dt, dC_dt, dP_dt]

# Solve the system numerically
initial_conditions = [A0, B0, C0, P0]  # Reordered for consistency
solution = odeint(reaction_rates, initial_conditions, t, args=(k1, k2))

# Extract the concentrations from the solution array
A_conc = solution[:, 0]
B_conc = solution[:, 1]
C_conc = solution[:, 2]
P_conc = solution[:, 3]  # Now C(t) is properly changing over time

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(t, A_conc, label=r'$A(t)$', linestyle='solid')
plt.plot(t, B_conc, label=r'$B(t)$', linestyle='dashdot')
plt.plot(t, C_conc, label=r'$C(t)$', linestyle='dotted', color="green")
plt.plot(t, P_conc, label=r'$P(t)$', linestyle='dashed', color="red")

# Labels and title
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.title('Kinetics of a Complex Successive Reaction')
plt.legend()
plt.grid(True)

# Show the plot
plt.show()