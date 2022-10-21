from ufl import (Coefficient, Constant, Measure, TestFunction, TrialFunction,
                 derivative, dot, grad)

from problem_specification import (T_fluid, T_remote, Vt, boltz, dom_map,
                                   emissivity, extrad_faces,
                                   flux, flux_faces, htc, mat_data, zone_faces)

# Get mesh through function space
mesh = Vt.ufl_domain()

# Define coefficients
T = Coefficient(Vt)
T0 = Coefficient(Vt)

v = TestFunction(Vt)
dT = TrialFunction(Vt)

# Define temporal variables
dt = Constant(mesh)

# Define integration measures
dx = Measure("dx", domain=mesh, metadata={"quadrature_degree": 4})
ds = Measure("ds", domain=mesh, metadata={"quadrature_degree": 4})

# Define variational form
F = 0
for dom in range(1, 6):
    material = mat_data[dom_map[dom]]
    kappa = material["kappa"]
    rho = material["rho"]
    cp = material["Cp"]
    # Add volume integrals
    F += rho*cp*(T - T0) * v * dx(dom) + dt * \
        dot(kappa * grad(T), grad(v)) * dx(dom)

# Add surface integrals
for ents in flux_faces:
    F += flux * v * ds(ents)
for ents in zone_faces:
    F += dt * htc * (T - T_fluid) * v * ds(ents)
for ents in extrad_faces:
    F += dt * boltz * emissivity * (T**4 - T_remote**4) * v * ds(ents)

J = derivative(F, T, dT)

forms = [F, J]
