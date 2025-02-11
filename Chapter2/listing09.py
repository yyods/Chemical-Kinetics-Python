from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# Define the system of differential equations for the parallel reaction.
# This function models a second-order parallel reaction:
#     CH3X + H2O  ->[k1] CH3OH   (Pathway 1)
#     CH3X + OH^-  ->[k2] CH3OH   (Pathway 2)
def parallel_reaction_kinetics(y, t, k1, k2, C_A0, C_B0):
    """
    Parameters:
        y    : list  -> [x], where x is the extent of reaction at time t.
        t    : float -> Time.
        k1   : float -> Rate constant for the hydrolysis pathway (with H2O).
        k2   : float -> Rate constant for the reaction with OH^-.
        C_A0 : float -> Initial concentration of CH3X.
        C_B0 : float -> Initial concentration of OH^-.
        
    Returns:
        A list containing the derivative dx/dt.
    """
    x = y[0]  # x(t) is the reacted fraction of CH3X.
    
    # Rate equation for the extent of reaction, considering both pathways.
    dx_dt = k1 * (C_A0 - x) + k2 * (C_A0 - x) * (C_B0 - x)
    
    return [dx_dt]

# Initial Conditions
C_A0 = 1.0  # Initial concentration of CH3X
C_B0 = 0.5  # Initial concentration of OH^-
x0   = [0]  # Initially, no reaction has occurred

# Rate Constants
k1 = 0.3  # Rate constant for hydrolysis with water
k2 = 0.1  # Rate constant for reaction with OH^-

# Time Grid for the Simulation
time = np.linspace(0, 10, 100)  # Simulate from t = 0 to t = 10 seconds

# Solve the ODE System Numerically
solution = odeint(parallel_reaction_kinetics, x0, time, args=(k1, k2, C_A0, C_B0))

# Extract the Concentration Profiles
x_t   = solution[:, 0]           # Extent of reaction
C_A_t = C_A0 - x_t               # Remaining CH3X concentration
# The following approximations assign portions of the reacted CH3X
# to the two pathways based on the relative rate constants.
C_B_t = (k1 / (k1 + k2 * C_B0)) * x_t   # CH3OH from the water pathway
C_C_t = (k2 * C_A0 / (k1 + k2 * C_B0)) * x_t  # CH3OH from the OH^- pathway

# Plot the Results
plt.figure(figsize=(8, 6))
plt.plot(time, C_A_t, label=r'$C_A(t)$ - CH3X', linestyle='solid', color='blue')
plt.plot(time, C_B_t, label=r'$C_B(t)$ - CH3OH (Water Pathway)', linestyle='dotted', color='green')
plt.plot(time, C_C_t, label=r'$C_C(t)$ - CH3OH (OH$^-$ Pathway)', linestyle='dashed', color='red')

# Graph Formatting
plt.xlabel('Time (t) [s]')
plt.ylabel('Concentration [M]')
plt.title('Second-Order Parallel Reaction Kinetics')
plt.legend()
plt.grid(True)
plt.show()