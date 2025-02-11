import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define symbols
t, gamma, lambda_, delta = sp.symbols('t gamma lambda delta', real=True, positive=True)
a0, b0, c0, d0, k1, k2 = sp.symbols('a0 b0 c0 d0 k1 k2', real=True, positive=True)

# Define x(t) as the extent of reaction
x = sp.Function('x')(t)

# Define the differential equation
ode = sp.Eq(sp.diff(x, t), -gamma * x**2 + lambda_ * x + delta)

# Solve the ODE with x(0)=0
sol = sp.dsolve(ode, x, ics={x.subs(t, 0): 0})
simplified_sol = sp.simplify(sol.rhs)

# Convert the solution to a numerical function
x_func = sp.lambdify((t, gamma, lambda_, delta), simplified_sol, 'numpy')

# Define parameter values
A0_val, B0_val, C0_val, D0_val = 0.06, 0.05, 0.04, 0.03
k1_val, k2_val = 0.5, 0.25

# Calculate gamma, lambda, and delta
gamma_val = k1_val + k2_val
lambda_val = k1_val * (A0_val + B0_val) - k2_val * (C0_val + D0_val)
delta_val = k1_val * A0_val * B0_val - k2_val * C0_val * D0_val

# Define time array
t_vals = np.linspace(0, 60, 300)
x_vals = x_func(t_vals, gamma_val, lambda_val, delta_val)

# Define numerical ODE model for the second-order reaction
def second_order_reversible(y, t, k1, k2, A0, B0, C0, D0):
    x = y[0]
    dx_dt = k1 * (A0 - x) * (B0 - x) - k2 * (C0 + x) * (D0 + x)
    return [dx_dt]

# Solve numerically with odeint
y0 = [0]  # x(0)=0
sol_x_numeric = odeint(second_order_reversible, y0, t_vals, args=(k1_val, k2_val, A0_val, B0_val, C0_val, D0_val))

# Compute concentrations
C_A_numeric = A0_val - sol_x_numeric[:, 0]
C_B_numeric = B0_val - sol_x_numeric[:, 0]
C_C_numeric = C0_val + sol_x_numeric[:, 0]
C_D_numeric = D0_val + sol_x_numeric[:, 0]

# Plot the concentration profiles
plt.figure(figsize=(8, 6))
plt.plot(t_vals, C_A_numeric, label=r'$A_0 - x(t)$', linestyle='solid', color='black')
plt.plot(t_vals, C_B_numeric, label=r'$B_0 - x(t)$', linestyle='dotted', color='black')
plt.plot(t_vals, C_C_numeric, label=r'$C_0 + x(t)$', linestyle='dashed', color='black')
plt.plot(t_vals, C_D_numeric, label=r'$D_0 + x(t)$', linestyle='dashdot', color='black')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.title('Second-Order Reversible Reaction Kinetics')
plt.legend()
plt.grid(True)
plt.show()