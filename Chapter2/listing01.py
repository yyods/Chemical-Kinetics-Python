from sympy import symbols, Eq, solve, integrate, pprint, latex, sqrt

# Define symbols for time (t), rate constant (k), and reaction order (n)
t, k, n = symbols('t k n')

# Define concentration variables: C_A (at time t) and initial concentration C_A0
C_A = symbols('C_A')
C_A0 = symbols('C_A0')

# Dictionary to store solutions for different reaction orders
solutions = {}

# Loop over different reaction orders: 0, 1, 2, and 3
for order in [0, 1, 2, 3]:
    # Separate variables: ∫ C_A^(-order) dC_A = -k ∫ dt
    lhs_integrated = integrate(C_A**(-order), C_A)
    rhs_integrated = integrate(-k, t)
    
    # Form the general integrated solution with an integration constant C1
    general_solution = Eq(lhs_integrated, rhs_integrated + symbols('C1'))
    
    # Use the initial condition C_A(0) = C_A0 to solve for C1
    C1_value = solve(general_solution.subs(t, 0).subs(C_A, C_A0), symbols('C1'))
    if C1_value:
        C1_value = C1_value[0]  # Extract the solution
        solution = solve(general_solution.subs(symbols('C1'), C1_value), C_A)
        if solution:
            solutions[order] = Eq(C_A, solution[0])  # Store the solution

# Display the solutions for each reaction order
for order, solution in solutions.items():
    print(f"Solution for reaction order n = {order}:")
    pprint(solution)
    print("LaTeX format:")
    print(latex(solution))
    print("\n" + "-"*50 + "\n")