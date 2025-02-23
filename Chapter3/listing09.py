import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Define symbolic variables
C_ES, C_E, C_E0, C_S, k1, k_neg1, k2 = sp.symbols('C_ES C_E C_E0 C_S k1 k_neg1 k2', positive=True, real=True)
KM = sp.symbols('KM', positive=True, real=True)

# Step 1: Equilibrium equation for ES complex formation
equilibrium_eq = sp.Eq(k1 * C_S * (C_E0 - C_ES), k_neg1 * C_ES)

# Step 2: Solve for C_ES
C_ES_solution = sp.solve(equilibrium_eq, C_ES)[0]

# Step 3: Define the rate equation
r_P_expr = k2 * C_ES_solution
sp.pprint(r_P_expr)

# Step 4: Express in terms of the Michaelis constant
r_P_MM = k2 * C_E0 * C_S / (C_S + k_neg1/k1)
r_P_final = sp.simplify(r_P_MM.subs(k_neg1/k1, KM))

print("The Michaelis–Menten rate expression is:")
sp.pretty_print(r_P_final)

# -----------------------------------------
# Numerical Evaluation and Plotting
# -----------------------------------------

rP_func = sp.lambdify((C_S, C_E0, k2, KM), r_P_final, 'numpy')

C_E0_val = 1.214e-7   # Total enzyme concentration (M)
k2_val   = 10         # Rate constant k2 (1/s)
KM_val   = 4.815e-4   # Michaelis constant (M)

C_S_values = np.linspace(1e-5, 0.015, 100)
rP_values = rP_func(C_S_values, C_E0_val, k2_val, KM_val)

# Plot 1: Reaction Rate vs. Substrate Concentration
plt.figure(figsize=(6, 4))
plt.plot(C_S_values, rP_values, 'k-', linewidth=2)
plt.axhline(y=k2_val * C_E0_val, linestyle='--', color='gray', label=r'$r_{max}$')
plt.xlabel(r'Substrate Concentration, $C_S$ (M)')
plt.ylabel(r'Reaction Rate, $r_P$ (M/s)')
plt.title('Reaction Rate vs. Substrate Concentration')
plt.legend()
plt.grid(True)
plt.show()

# Plot 2: Lineweaver–Burk Plot
inv_rP = 1 / rP_values    
inv_CS = 1 / C_S_values    
plt.figure(figsize=(6, 4))
plt.plot(inv_CS, inv_rP, 'k-', linewidth=2)
plt.xlabel(r'$1/C_S$ (1/M)')
plt.ylabel(r'$1/r_P$ (s/M)')
plt.title('Lineweaver–Burk Plot')
plt.grid(True)
plt.show()