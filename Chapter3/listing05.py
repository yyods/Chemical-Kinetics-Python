import sympy as sp

# 1) Define symbols for rate constants and concentrations.
#    Notation:
#      C1 : [C2H6]
#      C2 : [CH3*]
#      C3 : [CH4]
#      C4 : [C2H5*]
#      C5 : [H*]
#      C6 : [C2H4]
#      C7 : [H2]
#      C8 : [C4H10]
k1, k2, k3, k4, k5, k6 = sp.symbols('k1 k2 k3 k4 k5 k6', positive=True)
C1, C2, C3, C4, C5, C6, C7, C8 = sp.symbols('C1 C2 C3 C4 C5 C6 C7 C8', positive=True)

# 2) Define the reaction rates.
# Reaction scheme:
#  R1: C2H6 -> 2 CH3*
#  R2: CH3* + C2H6 -> CH4 + C2H5*
#  R3: C2H5* -> C2H4 + H*
#  R4: C2H6 + H* -> C2H5* + H2
#  R5: 2 C2H5* -> C4H10
#  R6: 2 C2H5* -> C2H4 + C2H6
r1 = k1 * C1                         # Reaction 1
r2 = k2 * C2 * C1                    # Reaction 2
r3 = k3 * C4                         # Reaction 3
r4 = k4 * C1 * C5                    # Reaction 4
r5 = k5 * C4**2                      # Reaction 5
r6 = k6 * C4**2                      # Reaction 6
r_vec = sp.Matrix([r1, r2, r3, r4, r5, r6])

# 3) Define the stoichiometric matrix.
# Rows correspond to the species:
#   0: C2H6, 1: CH3*, 2: CH4, 3: C2H5*, 4: H*, 5: C2H4, 6: H2, 7: C4H10
alpha = sp.Matrix([
    [-1, -1,  0, -1,  0, +1],   # C2H6
    [ 2, -1,  0,  0,  0,  0],   # CH3*
    [ 0, +1,  0,  0,  0,  0],   # CH4
    [ 0, +1, -1, +1, -2, -2],   # C2H5*
    [ 0,  0, +1, -1,  0,  0],   # H*
    [ 0,  0, +1,  0,  0, +1],   # C2H4
    [ 0,  0,  0, +1,  0,  0],   # H2
    [ 0,  0,  0,  0, +1,  0]    # C4H10
])

# 4) Compute the net production rates for each species.
net_rates = alpha * r_vec
# The entries of net_rates correspond to:
#   d[C2H6]/dt, d[CH3*]/dt, d[CH4]/dt, d[C2H5*]/dt, d[H*]/dt, d[C2H4]/dt, d[H2]/dt, d[C4H10]/dt

# 5) Apply the steady-state approximation for the radicals:
#    Set net rates for CH3* (row 1), C2H5* (row 3), and H* (row 4) equal to zero.
eq_CH3  = sp.Eq(net_rates[1], 0)  # CH3*
eq_C2H5 = sp.Eq(net_rates[3], 0)  # C2H5*
eq_H    = sp.Eq(net_rates[4], 0)  # H*
sol_radicals = sp.solve((eq_CH3, eq_C2H5, eq_H), (C2, C4, C5), dict=True)

print("Steady-state solutions for the radical species:")
sp.pprint(sol_radicals)
print("\n")
# Expected solutions:
#   [CH3*] = 2*k1/k2,
#   [C2H5*] = sqrt(k1 * C1/(k5 + k6)),
#   [H*]   = (k3/k4) * sqrt(k1/(C1*(k5+k6))).

# 6) Substitute the steady-state radical expressions into the net production rates
net_rates_ss = sp.simplify(net_rates.subs(sol_radicals[0]))

# 7) Display the differential equations for all species.
species = {
    0: r"[C2H6]",
    1: r"[CH3*]",
    2: r"[CH4]",
    3: r"[C2H5*]",
    4: r"[H*]",
    5: r"[C2H4]",
    6: r"[H2]",
    7: r"[C4H10]"
}

print("Differential equations after applying the steady-state approximation:\n")
for i in range(net_rates_ss.shape[0]):
    print(f"d{species[i]}/dt =")
    sp.pprint(net_rates_ss[i])
    print()