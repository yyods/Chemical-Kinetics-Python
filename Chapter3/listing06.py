import sympy as sp

# Define symbols and the function for ethane concentration C(t)
t, C0, alpha, beta = sp.symbols('t C0 alpha beta', positive=True)
C = sp.Function('C')(t)

# Original ODE: dC/dt = -alpha * C - beta * sqrt(C)
ode = sp.Eq(sp.diff(C, t), -alpha * C - beta * sp.sqrt(C))

# Substitute: let u(t) = sqrt(C(t)), so that C(t) = u(t)**2.
# Then, dC/dt = 2*u(t)*du/dt.
# Substituting into the ODE:
#   2*u*du/dt = -alpha*u**2 - beta*u.
# Dividing by u (assuming u>0) yields:
#   du/dt = - (alpha/2)*u - beta/2.
u = sp.Function('u')(t)
ode_u = sp.Eq(sp.diff(u, t), -alpha/2 * u - beta/2)

# Solve the ODE for u(t) with the initial condition: u(0) = sqrt(C0)
solution_u = sp.dsolve(ode_u, u, ics={u.subs(t, 0): sp.sqrt(C0)})

# The solution for C(t) is then [u(t)]**2.
solution_C = sp.simplify(solution_u.rhs**2)

print("Ethane concentration [C2H6](t) in terms of alpha and beta:")
sp.pprint(solution_C)