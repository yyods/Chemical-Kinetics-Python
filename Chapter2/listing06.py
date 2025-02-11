from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# Define the model for the basic consecutive reaction A -> B -> C
def model(y, t, k1, k2):
    A, B, C = y
    dA_dt = -k1 * A
    dB_dt = k1 * A - k2 * B
    dC_dt = k2 * B
    return [dA_dt, dB_dt, dC_dt]

# Initial conditions and parameters
A0, B0, C0 = 1.0, 0.0, 0.0
k1, k2 = 1.0/2.0, 1.0/3.0
time = np.linspace(0, 10, 100)

# Compute t_max for the intermediate B
t_max = np.log(k2 / k1) / (k2 - k1)

# Solve the ODE system
y0 = [A0, B0, C0]
sol = odeint(model, y0, time, args=(k1, k2))

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(time, sol[:, 0], label='[A]', linestyle='solid')
plt.plot(time, sol[:, 1], label='[B]', linestyle='dotted')
plt.plot(time, sol[:, 2], label='[C]', linestyle='dashed')
plt.axvline(x=t_max, linestyle="--", color="red", label=r"$t_{\max}$")
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.title('Consecutive Reaction Kinetics')
plt.legend()
plt.grid(True)
plt.show()