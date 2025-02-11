from sympy import symbols, Function, Eq, dsolve, pprint, exp, solve
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define symbols
t, k1, k2, A0, B0 = symbols('t k1 k2 A0 B0')
x = Function('x')(t)

# Define the differential equation
dx_dt = Eq(x.diff(t), k1 * (A0 - x) - k2 * (B0 + x))

# Solve the differential equation
sol_x = dsolve(dx_dt, x)
C1 = symbols('C1')

# Determine the integration constant using x(0)=0
C1_value = solve(sol_x.rhs.subs(t, 0) - 0, C1)[0]

# Substitute back to obtain x(t)
x_solution = sol_x.rhs.subs(C1, C1_value)

# Express concentrations in terms of x(t)
C_A = A0 - x_solution
C_B = B0 + x_solution

print("\nExplicit Solution for x(t):")
pprint(x_solution)
print("\nSolution for C_A(t):")
pprint(C_A)
print("\nSolution for C_B(t):")
pprint(C_B)

# Convert symbolic expressions to numerical functions
from sympy import lambdify
x_func = lambdify((t, k1, k2, A0, B0), x_solution, 'numpy')

# Numerical ODE model
def reversible_reaction(y, t, k1, k2, A0, B0):
    x = y[0]
    dx_dt = k1 * (A0 - x) - k2 * (B0 + x)
    return [dx_dt]

# Define parameters and initial conditions
A0_val, B0_val = 1.0, 0.2
k1_val, k2_val = 0.45, 0.12
t_vals = np.linspace(0, 10, 100)
y0 = [0]  # x(0)=0

# Solve numerically using odeint
sol_x_numeric = odeint(reversible_reaction, y0, t_vals, args=(k1_val, k2_val, A0_val, B0_val))

# Compute concentrations from numerical x(t)
C_A_numeric = A0_val - sol_x_numeric[:, 0]
C_B_numeric = B0_val + sol_x_numeric[:, 0]

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(t_vals, C_A_numeric, label=r'$C_A(t)$', linestyle='solid')
plt.plot(t_vals, C_B_numeric, label=r'$C_B(t)$', linestyle='dotted')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.title('Reversible Reaction Kinetics')
plt.legend()
plt.grid(True)
plt.show()