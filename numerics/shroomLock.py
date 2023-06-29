import dolfinx as dlfn
import meshio
from petsc4py import PETSc
import ufl
import os
import subprocess
from mpi4py import MPI
import pyvista
import numpy as np

# choose geometry:
geo_file = "meshes/Full-Lamelle.geo"

# create mesh using gmsh
# TODO see tutorial by dokken on how to use the gmsh API
msh_file = geo_file.replace(".geo", ".msh")
assert msh_file.endswith(".msh")
subprocess.run(["gmsh", "-v", "0", "-3", geo_file], stdout=subprocess.DEVNULL, check=True)

# read msh file
assert os.path.exists(msh_file)
mesh, cell_markers, facet_markers = dlfn.io.gmshio.read_from_msh(msh_file, MPI.COMM_WORLD, gdim=3)

# create function space
V = dlfn.fem.VectorFunctionSpace(mesh, ("CG", 2))
u, v = dlfn.fem.Function(V), ufl.TestFunction(V)
ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_markers)

## Define boundary conditions
# homogeneous Dirichlet:
u_D0 = np.array((0,) * mesh.geometry.dim, dtype=PETSc.ScalarType)
screw_dofs = dlfn.fem.locate_dofs_topological(V, facet_markers.dim, facet_markers.find(1))

# inhomogeneous Dirichlet:
#u_D1 = np.array((0., 0., 4.), dtype=PETSc.ScalarType)
u_Dz1 = PETSc.ScalarType(0.004)
inner_dofs = dlfn.fem.locate_dofs_topological(V.sub(2), facet_markers.dim, facet_markers.find(2))

# radial displacement
r_lamelle = 2.1 * 1e-3 # m #TODO nonlinear does not converge unless r_lamelle = 2.6mm
r_pilz = 3      * 1e-3 # m
radial_displacement = - (r_lamelle - r_pilz)

class myExpression:
    def __init__(self, r, u_r=0):
        self.r = r
        self.u_r = u_r
    def eval(self, x):
        return (self.u_r * x[0] / self.r,
                self.u_r *x[1] / self.r, 
                np.zeros_like(x[2]))
    def eval_x(self, x):
        return self.u_r * x[0] / self.r
    def eval_y(self, x):
        return  self.u_r *x[1] / self.r
    
u_Dr = myExpression(r_lamelle)
u_Dr.u_r = radial_displacement  
u_Dr_fun = dlfn.fem.Function(V)
#Vx = V.sub(0).collapse
#u_Dr_funx = dlfn.fem.Function(Vx)
#Vy = V.sub(1).collapse
#u_Dr_funy = dlfn.fem.Function(Vy)
u_Dr_fun.interpolate(u_Dr.eval)
#u_Dr_funx.interpolate(u_Dr.evalx)
#u_Dr_funy.interpolate(u_Dr.evaly)
u_Dr_funx, u_Dr_funy, u_Dr_funz = u_Dr_fun.split()

# FLäche für radialen Druck
inner_radial_dofs = dlfn.fem.locate_dofs_topological(V, facet_markers.dim, facet_markers.find(2))
inner_radial_dofsx = dlfn.fem.locate_dofs_topological(V.sub(0), facet_markers.dim, facet_markers.find(2))
inner_radial_dofsy = dlfn.fem.locate_dofs_topological(V.sub(1), facet_markers.dim, facet_markers.find(2))

# Fläche für flächigen Druck
r = r_pilz
def in_circle_bottom(x):
    return np.isclose(x[2], 0) * (np.sqrt(x[0] ** 2 + x[1] ** 2) < r)
def in_circle_top(x):
    return np.isclose(x[2], 1.2e-3) * (np.sqrt(x[0] ** 2 + x[1] ** 2) < r)

circle_facets_bottom = dlfn.mesh.locate_entities_boundary(mesh, mesh.topology.dim-1, in_circle_bottom)
circle_dofs_x_bottom = dlfn.fem.locate_dofs_topological(V.sub(0), mesh.topology.dim-1, circle_facets_bottom)
circle_dofs_y_bottom = dlfn.fem.locate_dofs_topological(V.sub(1), mesh.topology.dim-1, circle_facets_bottom)
circle_dofs_z_bottom = dlfn.fem.locate_dofs_topological(V.sub(2), mesh.topology.dim-1, circle_facets_bottom)

circle_facets_top = dlfn.mesh.locate_entities_boundary(mesh, mesh.topology.dim-1, in_circle_top)
circle_dofs_x_top = dlfn.fem.locate_dofs_topological(V.sub(0), mesh.topology.dim-1, circle_facets_top)
circle_dofs_y_top = dlfn.fem.locate_dofs_topological(V.sub(1), mesh.topology.dim-1, circle_facets_top)
circle_dofs_z_top = dlfn.fem.locate_dofs_topological(V.sub(2), mesh.topology.dim-1, circle_facets_top)


circle_dofs_bottom = dlfn.fem.locate_dofs_geometrical(V, in_circle_bottom)
circle_dofs_top = dlfn.fem.locate_dofs_geometrical(V, in_circle_top)


bcs = [dlfn.fem.dirichletbc(u_D0, screw_dofs, V),
       #dlfn.fem.dirichletbc(u_Dr_fun, inner_radial_dofs), # as u_Dr_fun is already interpolated, the dirichletbc does not need nor want any Info about V
       # Druck gegen Wand:
       #dlfn.fem.dirichletbc(u_Dr_funx, inner_radial_dofsx),
       #dlfn.fem.dirichletbc(u_Dr_funy, inner_radial_dofsy),
       # Druck auf Oberfläche:
       dlfn.fem.dirichletbc(u_Dr_funx, circle_dofs_x_top),
       dlfn.fem.dirichletbc(u_Dr_funy, circle_dofs_y_top),
       
]
       #dlfn.fem.dirichletbc(u_Dz1, circle_dofs_x, V.sub(2))]
       #dlfn.fem.dirichletbc(u_Dz1, inner_dofs, V.sub(2))]
       # dlfn.fem.dirichletbc(u_D1, inner_dofs, V)]

# Spatial dimension
d = len(u)

# Identity tensor
I = ufl.variable(ufl.Identity(d))

# Deformation gradient
F = ufl.variable(I + ufl.grad(u))

# Right Cauchy-Green tensor
C = ufl.variable(F.T * F)

# Invariants of deformation tensors
Ic = ufl.variable(ufl.tr(C))
J  = ufl.variable(ufl.det(F))

E = PETSc.ScalarType(40.0e6)
nu = PETSc.ScalarType(0.45)
mu = dlfn.fem.Constant(mesh, E/(2*(1 + nu)))
lmbda = dlfn.fem.Constant(mesh, E*nu/((1 + nu)*(1 - 2*nu)))

# Stored strain energy density (compressible neo-Hookean model)
psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J))**2
# Stress
# Hyper-elasticity
P = ufl.diff(psi, F)
# Linear-elasticity
#P = 2.0 * mu * ufl.sym(ufl.grad(u)) + lmbda * ufl.tr(ufl.sym(ufl.grad(u))) * I

metadata = {"quadrature_degree": 4}
dx = ufl.Measure("dx", domain=mesh, metadata=metadata)
# Neumann boundary condition
circle_tag = dlfn.mesh.meshtags(mesh, mesh.topology.dim - 1, circle_facets_bottom, 33)
ds = ufl.Measure("ds", domain=mesh, subdomain_data=circle_tag)
one = dlfn.fem.Constant(mesh, PETSc.ScalarType(1.))
f = dlfn.fem.form(one*ds(33))
A = dlfn.fem.assemble_scalar(f)
print(f"Contact surface of force: A = {A * 1e6} mm^2")
Force = 6.5 # N
traction = Force / A # N / m^2
t = dlfn.fem.Constant(mesh, PETSc.ScalarType((0, 0, -traction)))

F = ufl.inner(ufl.grad(v), P) * dx - ufl.inner(t, v) * ds(33)

problem = dlfn.fem.petsc.NonlinearProblem(F, u, bcs)

solver = dlfn.nls.petsc.NewtonSolver(mesh.comm, problem)
solver.atol = 1e-8
solver.rtol = 1e-8

solver.convergence_criterion = "incremental"
xdmf = dlfn.io.XDMFFile(MPI.COMM_WORLD, "Thilo--iterative-results_Lamelle.xdmf", "w")
xdmf.write_mesh(mesh)

# # TODO do we need this?!##################################
# V_coord = dlfn.fem.FunctionSpace(                       ##
#     mesh, mesh.ufl_domain().ufl_coordinate_element())   ##
# u_coord = dlfn.fem.Function(V_coord)                    ##
# ##########################################################

load_steps = 10
tVal = - traction / load_steps
u_DrVal = radial_displacement / load_steps


for n in range(1, load_steps + 1):
    print(f"------------------------------ Load step {n}: ------------------------------")
    # update displacement
    t.value[2] = n * tVal
    u_Dr.u_r = n * u_DrVal
    u_Dr_fun.interpolate(u_Dr.eval)
    u_Dr_funx, u_Dr_funy, u_Dr_funz = u_Dr_fun.split()
    num_its, converged = solver.solve(u)
    u.name = "Deflection"
    xdmf.write_function(u, n)
    u.x.scatter_forward() # TODO this has something to do with MPI stuff which I do not understand 
    assert(converged)
    print(f"Load step {n}, Iterations {num_its}, F = {np.abs(t.value[2] * A):.2f} N, u_r = {u_Dr.u_r * 1e3:.2f} mm, max(|u_z|) = {np.max(np.abs(u.x.array.reshape(-1 , len(u))[:,2])) * 1e3:.2f} mm\n")
    # # TODO do we need this?!######################################################################
    # u_coord.interpolate(u)                                                                      ##
    # mesh.geometry.x[:, :mesh.geometry.dim] += u_coord.x.array.reshape((-1, mesh.geometry.dim))  ##
    # one = dlfn.fem.Constant(mesh, PETSc.ScalarType(1.))                                         ##
    # f = dlfn.fem.form(one*ds(33))                                                               ##
    # A = dlfn.fem.assemble_scalar(f)                                                             ##
    # print(f"A = {A * 1e-6}  mm^2")                                                              ##
    # ##############################################################################################


## post processing # TODO I do not trust my implementation of the hyper elasticity
# Hyper-elasticity
# Deformation gradient
F = ufl.variable(I + ufl.grad(u))

# Right Cauchy-Green tensor
C = ufl.variable(F.T * F)

# Invariants of deformation tensors
Ic = ufl.variable(ufl.tr(C))
J  = ufl.variable(ufl.det(F))

# Stored strain energy density (compressible neo-Hookean model)
psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J))**2
# Stress
P = ufl.diff(psi, F)

# Linear-elasticity
# stress
#P = 2.0 * mu * ufl.sym(ufl.grad(u)) + lmbda * ufl.tr(ufl.sym(ufl.grad(u))) * I

# get von Mises
vM = ufl.sqrt(2./3. * ufl.inner(ufl.tr(P), ufl.tr(P)))
V_vonMises = dlfn.fem.FunctionSpace(mesh, ("DG", 0)) # TODO why DG0, not CG1?
vM_expr = dlfn.fem.Expression(vM, V_vonMises.element.interpolation_points())
vM_Fun = dlfn.fem.Function(V_vonMises)
vM_Fun.interpolate(vM_expr)

vM_Fun.name = "vonMises"
#with dlfn.io.XDMFFile(MPI.COMM_WORLD, "wild-Testing-results_Lamelle.xdmf", "w") as xdmf:
#xdmf.write_mesh(mesh)
#    xdmf.write_function(u, 0)
xdmf.write_function(vM_Fun, n)

# # plotting
# topology, cell_types, x = dlfn.plot.create_vtk_mesh(mesh, mesh.topology.dim)
# grid = pyvista.UnstructuredGrid(topology, cell_types, x)
# p = pyvista.Plotter(window_size=[800, 800])
# p.add_mesh(grid, show_edges=True)
# p.show()