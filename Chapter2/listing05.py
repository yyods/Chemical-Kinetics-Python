from sympy import symbols, Function, Eq, dsolve, pprint, exp, solve

# Define symbols and functions
t, k1, k2, C_A0 = symbols('t k1 k2 C_A0')
C_A = Function('C_A')(t)
C_B = Function('C_B')(t)
C_C = Function('C_C')(t)

# Differential equations for the basic mechanism
eq1 = Eq(C_A.diff(t), -k1 * C_A)
eq2 = Eq(C_B.diff(t), k1 * C_A - k2 * C_B)

# Solve for C_A(t)
sol_A = dsolve(eq1, C_A)
C_A_solution = sol_A.rhs.subs('C1', C_A0)

# Substitute C_A(t) into eq2 and solve for C_B(t)
eq2_subs = eq2.subs(C_A, C_A_solution)
sol_B = dsolve(eq2_subs, C_B)
C1_B = symbols('C1_B')
general_solution_B = sol_B.rhs.subs('C1', C1_B)
C1_B_value = solve(general_solution_B.subs(t, 0), C1_B)[0]
C_B_solution = general_solution_B.subs(C1_B, C1_B_value)

print("Solution for C_A(t):")
pprint(C_A_solution)
print("\nSolution for C_B(t):")
pprint(C_B_solution)