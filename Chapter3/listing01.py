import sympy as sp

# Define symbols
A, B, C, D, E = sp.symbols('A B C D E')
k1, k2, k3, k4, k5, k6 = sp.symbols('k1 k2 k3 k4 k5 k6')

# Define rate equations
r = sp.Matrix([
    k1*A**2, 
    k2*B**2*C, 
    k3*A*B, 
    k4*D**2, 
    k5*B*D, 
    k6*E
])

# Stoichiometric matrix
alpha = sp.Matrix([[-2, 2, 1, 0, 0],
                    [ 2, -2, -1, 0, 0],
                    [-1, -1, 0, 2, 0],
                    [ 1, 1, 0, -2, 0],
                    [ 0, -1, 0, -1, 1],
                    [ 0, 1, 0, 1, -1]])

# Compute concentration derivatives
dC_dt = alpha.T * r
sp.pprint(dC_dt)