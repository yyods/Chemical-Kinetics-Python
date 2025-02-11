import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# 1) SYMBOLIC DERIVATION: Solve for [A*] under the steady-state condition.
#    The steady-state condition is: d[A*]/dt = 0, i.e.,
#       k1*[A]^2 = (k2*[A] + k3)*[A*]
# ---------------------------------------------------------------

# Define symbolic variables
A_sym, Astar_sym = sp.symbols('A Astar', positive=True)
k1_sym, k2_sym, k3_sym = sp.symbols('k1 k2 k3', positive=True)

# Write the steady-state equation for A*
eq_ss = sp.Eq(k1_sym*A_sym**2, Astar_sym * (k2_sym*A_sym + k3_sym))
# Solve for [A*]
Astar_ss_expr = sp.solve(eq_ss, Astar_sym)[0]

print("Symbolic steady-state solution for [A*]:")
sp.pprint(Astar_ss_expr)

# ---------------------------------------------------------------
# 2) NUMERICAL INTEGRATION: Solve the full ODE system.
#
#    d[A]/dt    = -k1*A^2 + k2*A*A*
#    d[A*]/dt   =  k1*A^2 - k2*A*A* - k3*A*
#    d[P]/dt    =  k3*A*
# ---------------------------------------------------------------

# Define numerical values for the rate constants
k1_val = 1.0
k2_val = 10.0
k3_val = 1.0

# Initial conditions: [A]=1.0, [A*]=0.0, [P]=0.0
A0 = 1.0
Astar0 = 0.0
P0 = 0.0

def lindemann_odes(t, y):
    A, Astar, P = y
    dA_dt     = -k1_val*A**2 + k2_val*A*Astar
    dAstar_dt =  k1_val*A**2 - k2_val*A*Astar - k3_val*Astar
    dP_dt     =  k3_val*Astar
    return [dA_dt, dAstar_dt, dP_dt]

# Time span and points for evaluation
t_span = (0, 5)
t_eval = np.linspace(0, 5, 200)
sol = solve_ivp(lindemann_odes, t_span, [A0, Astar0, P0], t_eval=t_eval)
t_vals = sol.t
A_vals = sol.y[0]
Astar_vals = sol.y[1]
P_vals = sol.y[2]

# ---------------------------------------------------------------
# 3) EVALUATE THE STEADY-STATE EXPRESSION:
#    Convert the symbolic expression for [A*] into a function and
#    evaluate it using the numerical [A] from the ODE solution.
# ---------------------------------------------------------------

Astar_ss_func = sp.lambdify(
    (A_sym, k1_sym, k2_sym, k3_sym),
    Astar_ss_expr,
    'numpy'
)
Astar_ss_vals = Astar_ss_func(A_vals, k1_val, k2_val, k3_val)

# ---------------------------------------------------------------
# 4) PLOTTING THE RESULTS:
#    Compare the numerical solution for [A*] with the steady-state
#    approximation.
# ---------------------------------------------------------------
plt.figure(figsize=(8,6))
plt.plot(t_vals, A_vals,       'b-', label='[A](t) (numeric)')
plt.plot(t_vals, Astar_vals,   'r-', label='[A*](t) (numeric)')
plt.plot(t_vals, P_vals,       'g-', label='[P](t) (numeric)')
plt.plot(t_vals, Astar_ss_vals,'r--', label='[A*](t) steady-state')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.title('Lindemann Mechanism: Comparison with Steady-State Approximation')
plt.legend()
plt.grid(True)
plt.show()