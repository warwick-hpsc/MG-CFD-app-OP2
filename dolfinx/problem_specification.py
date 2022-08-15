# Copyright (C) 2022 Jørgen S. Dokken
#
# SPDX-License-Identifier:    MIT
#
from ufl import FiniteElement, FunctionSpace, Mesh, VectorElement, tetrahedron

# Define domain (mesh)
cell = tetrahedron
c_el = VectorElement("Lagrange", cell, 1, dim=3)
mesh = Mesh(c_el)

# Define thermal function space
t_element = FiniteElement("Lagrange", cell, 1)
Vt = FunctionSpace(mesh, t_element)

# Define structural function space
s_element = VectorElement("Lagrange", cell, 1)
Vs = FunctionSpace(mesh, s_element)


#model = read_model_json("../crescendo-mesh/CRESCENDO_ENGINE_LTMS_order1_cpp.json", {"T": T, "time": time}, r)
# Cell information
mat_data = {"TNM": {"rho": 0.443000E-08, "nu": 0.33, "E": 116000.0, "alpha": 0.852000E-05, "kappa": 6.77000, "Cp": 0.545000E+09},
            "NKL": {"rho": 0.822000E-08, "nu": 0.32, "E": 207000.0, "alpha": 0.129000E-04, "kappa": 10.9000, "Cp": 0.424000E+09},
            "STL": {"rho": 0.777000E-08, "nu": 0.28, "E": 212000.0, "alpha": 0.979000E-05, "kappa": 20.9000, "Cp": 0.449000E+09}
            }
dom_map = {1: "NKL", 2: "TNM", 3: "NKL", 4: "NKL", 5: "STL"}

# Face information
flux_faces = [512]
zone_faces = [517]
extrad_faces = [521, 526]
pres_faces = [521, 526]
tl_faces = []


# Variables used in Transient-Thermal simulation
htc = 0.45
flux = 10.0
T_remote = 293.0
T_fluid = 393.0
emissivity = 0.05
boltz = 5.6703E-18

# Variables used in Structural simulation
T_ref = 293.15
PRES = 100.0
TL = 1000.0
