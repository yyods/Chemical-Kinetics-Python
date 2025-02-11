import sympy as sp

# Define symbols
t = sp.Symbol('t')
C0 = sp.Matrix([sp.Symbol('C0_1'), 0, 0, 0])
k1, k2, k3, k4 = sp.symbols('k1 k2 k3 k4')

# Define the rate constant matrix
K = sp.Matrix([[-(k1 + k2), 0, 0, 0],
               [k1, 0, 0, 0],
               [k2, 0, -k3, k4],
               [0, 0, k3, -k4]])

# Compute eigenvalues and eigenvectors
eigen_data = K.eigenvects()
eigenvalues = [ev[0] for ev in eigen_data]
eigenvectors = []
for ev in eigen_data:
    for vec in ev[2]:  # Include all eigenvectors for degenerate cases
        eigenvectors.append(vec)

# Ensure the eigenvector matrix is square
if len(eigenvectors) == K.shape[0]:
    X = sp.Matrix.hstack(*eigenvectors)  # Construct eigenvector matrix
    e_Lambda_t = sp.diag(*[sp.exp(ev[0] * t) for ev in eigen_data for _ in ev[2]])
    
    # Compute matrix exponential using similarity transformation
    exp_Kt = X * e_Lambda_t * X.inv()

    # Compute concentration evolution
    C_t = exp_Kt * C0

    # Display final concentration evolution
    print("\nConcentration Evolution:")
    sp.pprint(C_t)
else:
    print("Error: Eigenvector matrix is not square. Direct diagonalization is not possible.")