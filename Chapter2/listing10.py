from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp  # Import sympy for symbolic computations

# Define the ODE for the self-catalyzed reaction
def self_catalyzed(y, t, k, CA0, CB0):
    x = y[0]  # Reacted fraction of acetone
    dx_dt = k * (CA0 - x) * (CB0 + x)
    return [dx_dt]

# Initial conditions
CA0, CB0 = 0.8, 0.001  # Initial concentrations of acetone and hydronium ions
x0 = [0]              # Initial reacted fraction
k = 0.2               # Rate constant

# Define the time points for integration
time = np.linspace(0, 100, 200)

# Solve the differential equation
solution = odeint(self_catalyzed, x0, time, args=(k, CA0, CB0))

# Extract the data
x_t = solution[:, 0]
CB_t = CB0 + x_t  # Concentration of hydronium ions (product B)
CA_t = CA0 - x_t  # Concentration of acetone
reaction_rate = k * CA_t * CB_t  # Reaction rate

# Symbolic computation of t_max using sympy
t_sym, CA0_sym, CB0_sym, k_sym = sp.symbols('t CA0 CB0 k')
t_max_expr = sp.log(CA0_sym / CB0_sym) / (k_sym * (CA0_sym + CB0_sym))  # Expression for t_max

# Pretty print the symbolic expression for t_max
print("\nSymbolic expression for t_max:")
sp.pprint(t_max_expr)

# Evaluate t_max numerically
t_max_value = t_max_expr.subs({CA0_sym: CA0, CB0_sym: CB0, k_sym: k}).evalf()
print(f"\nNumerical value of t_max = {t_max_value:.4f}")

# Create the figure and subplots
fig, axes = plt.subplots(2, 1, figsize=(8, 10))

# First plot: Concentration profiles
axes[0].plot(time, CA_t, label=r'$C_A(t)$', linestyle='solid')
axes[0].plot(time, CB_t, label=r'$C_B(t)$', linestyle='dashed')
axes[0].set_xlabel('Time (t)')
axes[0].set_ylabel('Concentration')
axes[0].set_title('Self-Catalyzed Reaction Kinetics')
axes[0].legend()
axes[0].grid(True)

# Second plot: Reaction rate vs. time
axes[1].plot(time, reaction_rate, label=r'$k C_A(t) C_B(t)$', linestyle='solid')
axes[1].axvline(x=float(t_max_value), color='red', linestyle='dashed', label=r'$t_{\max}$')
axes[1].set_xlabel('Time (t)')
axes[1].set_ylabel('Reaction Rate')
axes[1].set_title('Reaction Rate vs Time')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()