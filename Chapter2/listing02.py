import numpy as np
import matplotlib.pyplot as plt

# Define numerical parameters
k_value = 1      # Rate constant
C_A0_value = 1   # Initial concentration
time = np.linspace(0, 1, 100)  # Time points

# Define analytical functions for each reaction order
def C_A_0order(t):
    return np.maximum(C_A0_value - k_value * t, 0)

def C_A_1order(t):
    return C_A0_value * np.exp(-k_value * t)

def C_A_2order(t):
    return C_A0_value / (1 + k_value * C_A0_value * t)

def C_A_3order(t):
    return C_A0_value / np.sqrt(1 + 2 * k_value * C_A0_value**2 * t)

# Plot the curves
plt.figure(figsize=(8, 6))
plt.plot(time, C_A_0order(time), label='n = 0', linestyle='solid')
plt.plot(time, C_A_1order(time), label='n = 1', linestyle='dotted')
plt.plot(time, C_A_2order(time), label='n = 2', linestyle='dashed')
plt.plot(time, C_A_3order(time), label='n = 3', linestyle='dashdot')

plt.xlabel('Time (t)')
plt.ylabel('Concentration $C_A$')
plt.title('Kinetic Curves for Different Reaction Orders')
plt.legend()
plt.grid(True)
plt.show()