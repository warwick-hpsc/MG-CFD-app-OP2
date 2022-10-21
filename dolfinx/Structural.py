from ufl import (Coefficient, Constant, FacetNormal, Identity,
                 SpatialCoordinate, TestFunction, TrialFunction, derivative,
                 dot, ds, dx, grad, inner, sqrt, tr)

from problem_specification import (PRES, T_ref, Vs, Vt, dom_map, mat_data,
                                   pres_faces)

# Get mesh through function space
mesh = Vs.ufl_domain()

# Define spatial variables
n = FacetNormal(mesh)
x = SpatialCoordinate(mesh)
r = sqrt(x[1]**2 + x[2]**2)

# Define coefficients
T = Coefficient(Vt)
u = Coefficient(Vs)

du = TrialFunction(Vs)
v = TestFunction(Vs)

# Define temporal variables
time = Constant(mesh)

# Define helper functions


def epsilon(u):
    return 0.5*(grad(u) + grad(u).T)


def epsilon_t(t, alpha):
    return alpha*(t - T_ref)*Identity(3)


def sigma(u, t, alpha):
    eps = epsilon(u) - epsilon_t(t, alpha)
    return 2.0*mu*eps + lmbda*tr(eps)*Identity(3)


# Define variational form
F = 0
for dom in range(1, 6):
    material = mat_data[dom_map[dom]]
    alpha = material["alpha"]
    E = material["E"]
    nu = material["nu"]

    mu = E/(2.0*(1.0 + nu))
    lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))

    F += inner(sigma(u, T, alpha), epsilon(v)) * dx(dom)


for ents in pres_faces:
    F -= dot(PRES*n, v) * ds(ents)

# One = Constant(mesh)
# for ents in tl_faces:
#     area = One*ds(ents)
#     F -= dot((TL/area)*n, v) * ds(ents)


J = derivative(F, u, du)

forms = [J, F]
