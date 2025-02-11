from sympy import symbols, Eq, solve, pprint

# Define symbols
k, A = symbols('k A')

# Define a simple rate equation
rate_eq = Eq(A - k*A, 1)
pprint(rate_eq)

# Solve for A
solution = solve(rate_eq, A)
pprint(solution)