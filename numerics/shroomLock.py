import dolfinx as dlfn
from petsc4py import PETSc
from petsc4py.PETSc import ScalarType
import ufl
import os
import subprocess
from mpi4py import MPI
import numpy as np

# choose geometry:
geo_file = "meshes/ShroomLock-V7-Lock-quarter.geo"
#geo_file = "meshes/ShroomLock-V7-Lock.geo"


# create mesh using gmsh
# TODO use the gmsh API, for now a local gmsh version is used with subprocess
msh_file = geo_file.replace(".geo", ".msh")
assert msh_file.endswith(".msh")
subprocess.run(["gmsh", "-v", "0", "-3", geo_file], stdout=subprocess.DEVNULL, check=True)

# read msh file
assert os.path.exists(msh_file)
mesh, cell_markers, facet_markers = dlfn.io.gmshio.read_from_msh(msh_file, MPI.COMM_WORLD, gdim=3)
normal = ufl.FacetNormal(mesh)

# create function space
V = dlfn.fem.VectorFunctionSpace(mesh, ("CG", 2))
u, v = dlfn.fem.Function(V), ufl.TestFunction(V)
du = ufl.TrialFunction(V)

# define values
u_D0 = dlfn.fem.Constant(mesh, PETSc.ScalarType((0, 0, 0)))
r_lamelle = 2.6 * 1e-3 # m
r_pilz = 3.4      * 1e-3 # m
radial_displacement_value = r_pilz - r_lamelle
thickness_lamelle = 1.5 * 1e-3 # m

class radial_displacement_class:
    def __init__(self, r, u_r=0):
        self.r = r
        self.u_r = u_r
    def eval(self, x):
        return (self.u_r * x[0] / self.r,
                self.u_r * x[1] / self.r, 
                np.zeros_like(x[2]))
radial_displacement = radial_displacement_class(r_lamelle)
radial_displacement.u_r = radial_displacement_value
radial_displacement_fun = dlfn.fem.Function(V)
radial_displacement_fun.interpolate(radial_displacement.eval)
radial_displacement_funx, radial_displacement_funy, radial_displacement_funz = radial_displacement_fun.split()

dim = len(u)                        # Spatial dimension
I = ufl.variable(ufl.Identity(dim)) # Identity tensor
F = ufl.variable(I + ufl.grad(u))   # Deformation gradient
C = ufl.variable(F.T * F)           # Right Cauchy-Green tensor
Ic = ufl.variable(ufl.tr(C))        # Invariants of deformation tensors
J  = ufl.variable(ufl.det(F))       # Invariants of deformation tensors

E = PETSc.ScalarType(40.0e6)
nu = PETSc.ScalarType(0.45)
mu = dlfn.fem.Constant(mesh, E/(2*(1 + nu)))
lmbda = dlfn.fem.Constant(mesh, E*nu/((1 + nu)*(1 - 2*nu)))


psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J))**2 # Stored strain energy density (compressible neo-Hookean model)
P = ufl.diff(psi, F)                                                      # Stress for compressible neo-Hookean model

#P = 2.0 * mu * ufl.sym(ufl.grad(u)) + lmbda * ufl.tr(ufl.sym(ufl.grad(u))) * I # Stress for Linear-elasticity

## Define boundaries
screw_dofs = dlfn.fem.locate_dofs_topological(V, facet_markers.dim, facet_markers.find(1))
inner_radial_dofs = dlfn.fem.locate_dofs_topological(V, facet_markers.dim, facet_markers.find(2))
inner_radial_dofsx = dlfn.fem.locate_dofs_topological(V.sub(0), facet_markers.dim, facet_markers.find(2))
inner_radial_dofsy = dlfn.fem.locate_dofs_topological(V.sub(1), facet_markers.dim, facet_markers.find(2))
inner_radial_dofsz = dlfn.fem.locate_dofs_topological(V.sub(2), facet_markers.dim, facet_markers.find(2))
sym_xx_dofs = dlfn.fem.locate_dofs_topological(V.sub(0), facet_markers.dim, facet_markers.find(101))
sym_xy_dofs = dlfn.fem.locate_dofs_topological(V.sub(1), facet_markers.dim, facet_markers.find(101))
sym_yx_dofs = dlfn.fem.locate_dofs_topological(V.sub(0), facet_markers.dim, facet_markers.find(102))
sym_yy_dofs = dlfn.fem.locate_dofs_topological(V.sub(1), facet_markers.dim, facet_markers.find(102))

# surface for contacts in top and bottom of lamelle
r = r_pilz
def in_circle_bottom(x):
    return np.isclose(x[2], 0) * (np.sqrt(x[0] ** 2 + x[1] ** 2) < r) # TODO doesnt work with Lip in v7
circle_facets_bottom = dlfn.mesh.locate_entities_boundary(mesh, mesh.topology.dim-1, in_circle_bottom)
circle_dofs_x_bottom = dlfn.fem.locate_dofs_topological(V.sub(0), mesh.topology.dim-1, circle_facets_bottom)
circle_dofs_y_bottom = dlfn.fem.locate_dofs_topological(V.sub(1), mesh.topology.dim-1, circle_facets_bottom)
circle_dofs_z_bottom = dlfn.fem.locate_dofs_topological(V.sub(2), mesh.topology.dim-1, circle_facets_bottom)

def in_circle_top(x):
    return np.isclose(x[2], thickness_lamelle) * (np.sqrt(x[0] ** 2 + x[1] ** 2) < r_pilz)
def in_circle_radial(x):
    return np.isclose(np.sqrt(x[0] ** 2 + x[1] ** 2), r_lamelle) * (x[2] >= 0) * (x[2] <= thickness_lamelle)
def on_top(x):
    return np.isclose(x[2], 6 * 1e-3)
def on_bottom(x):
    return np.isclose(x[2], 0)

circle_facets_top = dlfn.mesh.locate_entities_boundary(mesh, mesh.topology.dim-1, in_circle_top)

circle_facets_inner = dlfn.mesh.locate_entities_boundary(mesh, mesh.topology.dim-1, in_circle_radial)
top_facets = dlfn.mesh.locate_entities(mesh, mesh.topology.dim-1, on_top)
bottom_facets = dlfn.mesh.locate_entities(mesh, mesh.topology.dim-1, on_bottom)
marked_facets = np.hstack([circle_facets_top,
                           circle_facets_inner,
                           top_facets])

marked_values = np.hstack([np.full_like(circle_facets_top, 4),
                           np.full_like(circle_facets_inner, 2), 
                           np.full_like(top_facets, 66)])

sorted_facets = np.argsort(marked_facets)
facet_tag = dlfn.mesh.meshtags(mesh, mesh.topology.dim-1, marked_facets[sorted_facets], marked_values[sorted_facets])

top_dofs = dlfn.fem.locate_dofs_topological(V, mesh.topology.dim-1, top_facets)
bottom_dofs = dlfn.fem.locate_dofs_topological(V, mesh.topology.dim-1, bottom_facets)

circle_dofs_x_top = dlfn.fem.locate_dofs_topological(V.sub(0), mesh.topology.dim-1, circle_facets_top)
circle_dofs_y_top = dlfn.fem.locate_dofs_topological(V.sub(1), mesh.topology.dim-1, circle_facets_top)
circle_dofs_z_top = dlfn.fem.locate_dofs_topological(V.sub(2), mesh.topology.dim-1, circle_facets_top)
def ppos(x):
    return (x + abs(x)) / 2
def pneg(x):
    return (x - abs(x)) / 2
# obstacle (penalty approach for contact)
class obstacle_class:
    def __init__(self, r, u_r=0, z_displacement=0.0):
        self.r = r
        self.u_r = u_r
        self.z_displacement = z_displacement
        
    def eval(self, x):
        rxy = np.sqrt(x[0] ** 2 + x[1] ** 2)
        
        r1 = 2.4 * 1e-3
        r2 = 3.4  * 1e-3
        h1 = 1.54 * 1e-3
        h2 = 3.5  * 1e-3
        u_r = (r2 - r_lamelle) * self.z_displacement / (4 * 1e-3) # should ideally be dependent on u[2] solution component
        return (
                u_r * x[0] / rxy,
                u_r * x[1] / rxy,
                - self.z_displacement + ((h2 - h1) / (r2 - r1) * (rxy - r1))# * (rxy < r2) + 4 * (h2 - h1) * np.ones_like(rxy) * (rxy >= r2)
                )

obstacle = obstacle_class(r_lamelle)
obstacle_fun = dlfn.fem.Function(V)
obstacle_fun.interpolate(obstacle.eval)
obstacle_funx, obstacle_funy, obstacle_funz = obstacle_fun.split()

circle_tag = dlfn.mesh.meshtags(mesh, mesh.topology.dim - 1, circle_facets_top, 4)
#-------------------------------------------------------------------------------------#
mesh.topology.create_connectivity(mesh.topology.dim-1, mesh.topology.dim)
with dlfn.io.XDMFFile(mesh.comm, "results/facet_tags.xdmf", "w") as xdmf:
    xdmf.write_mesh(mesh)
    xdmf.write_meshtags(facet_tag)
    #xdmf.write_meshtags(facet_markers)
#-------------------------------------------------------------------------------------#

# define boundary conditions
bcs = [dlfn.fem.dirichletbc(u_D0, top_dofs, V),
       dlfn.fem.dirichletbc(ScalarType(0), sym_xx_dofs, V.sub(0)),
       dlfn.fem.dirichletbc(ScalarType(0), sym_xy_dofs, V.sub(1)),
       dlfn.fem.dirichletbc(ScalarType(0), sym_yx_dofs, V.sub(0)),
       dlfn.fem.dirichletbc(ScalarType(0), sym_yy_dofs, V.sub(1)),
       #dlfn.fem.dirichletbc(ScalarType(0), circle_dofs_z_top, V.sub(2)) 
    ]  

metadata = {"quadrature_degree": 4}
dx = ufl.Measure("dx", domain=mesh, metadata=metadata)
#ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_markers)
ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tag)
#ds1 = ufl.Measure("ds", domain=mesh, subdomain_data=circle_tag)

# weak form
pen = dlfn.fem.Constant(mesh, ScalarType(1e14))
A_screws = dlfn.fem.assemble_scalar(dlfn.fem.form(dlfn.fem.Constant(mesh, ScalarType(1)) * ds(1))) # ds(1) is not defined anymore
print(f"Top surface area: {dlfn.fem.assemble_scalar(dlfn.fem.form(dlfn.fem.Constant(mesh, ScalarType(1)) * ds(66)))}")
r1 = 2.4 * 1e-3
r2 = 3.4  * 1e-3
h1 = 1.54 * 1e-3
h2 = 3.5  * 1e-3
t = 1.5  * 1e-3

form = ufl.inner(ufl.grad(v), P) * dx 
# form += pen * (ufl.dot(v[0], u[0] - obstacle_funx)) * ds(41)
# form += pen * (ufl.dot(v[1], u[1] - obstacle_funy)) * ds(41)
#----------------------------------------------------------------------------------------------------------------------------------------------------
# holy, this is actually working!
x = ufl.SpatialCoordinate(mesh)
rxy = ufl.sqrt(x[0]**2 + x[1]**2)
# form += pen * (ufl.dot(v[0], ppos(-(0.5 * 1e-3) * x[0] / r_lamelle + u[0]))) * ds(4) # we are on the first quadrant (obstacle and u are negative)
# form += pen * (ufl.dot(v[1], pneg(-(0.5 * 1e-3) * x[1] / r_lamelle + u[1]))) * ds(4) # we are on the first quadrant (obstacle and u are negative)
# #----------------------------------------------------------------------------------------------------------------------------------------------------
# form += pen * (ufl.dot(v[0], ppos((r2 - r_lamelle) / (h2 - h1) * u[2] * x[0] / rxy + u[0]))) * ds(2) # we are on the first quadrant (obstacle and u are negative)
# form += pen * (ufl.dot(v[1], pneg((r2 - r_lamelle) / (h2 - h1) * u[2] * x[1] / rxy + u[1]))) * ds(2) # we are on the first quadrant (obstacle and u are negative)
# form += pen * (ufl.dot(v[0], ppos(- obstacle_funx + u[0]))) * ds(4) # we are on the first quadrant (obstacle and u are negative)
# form += pen * (ufl.dot(v[1], pneg(- obstacle_funy + u[1]))) * ds(4) # we are on the first quadrant (obstacle and u are negative)
form += pen * (ufl.dot(v[2], ppos(u[2] - obstacle_funz))) * ds(4) # this is the correct one
# form += pen * ufl.dot(v[0], u[0] - obstacle_funx) * ds(4)
# form += pen * ufl.dot(v[1], u[1] - obstacle_funy) * ds(4)

# Force approach ------------------------------------------------------------------------------------------------------------------------------------
# normal_traction = dlfn.fem.Constant(mesh, PETSc.ScalarType(0))
# tangential_traction = dlfn.fem.Constant(mesh, PETSc.ScalarType(0))
# form += ufl.dot(v[2], normal_traction) * ds(4)
# form += ufl.dot(v[0], - tangential_traction * x[0] / rxy) * ds(4)
# form += ufl.dot(v[1], - tangential_traction * x[1] / rxy) * ds(4)

J_form = ufl.derivative(form, u, du)
problem = dlfn.fem.petsc.NonlinearProblem(form, u, bcs, J=J_form)

solver = dlfn.nls.petsc.NewtonSolver(mesh.comm, problem)
solver.atol = 1e-8  
solver.rtol = 1e-8

solver.convergence_criterion = "incremental"
xdmf = dlfn.io.XDMFFile(MPI.COMM_WORLD, "results/quarter-results_simpleContact-nonlin.xdmf", "w")
xdmf.write_mesh(mesh)

## get Force on top
def get_force(u):
    F = ufl.variable(I + ufl.grad(u))
    C = ufl.variable(F.T * F)
    
    Ic = ufl.variable(ufl.tr(C))
    J  = ufl.variable(ufl.det(F))

    psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J))**2
    P = ufl.diff(psi, F)
    
    #P = 2.0 * mu * ufl.sym(ufl.grad(u)) + lmbda * ufl.tr(ufl.sym(ufl.grad(u))) * I
    
    n_top = dlfn.fem.Constant(mesh, np.array([0,0,1], dtype=np.float64))
    trac = ufl.inner(P, ufl.outer(n_top, n_top)) * ds(66)

    return dlfn.fem.assemble_scalar(dlfn.fem.form(trac))

## get von Mises
def get_vonMises(u):
    F = ufl.variable(I + ufl.grad(u))
    C = ufl.variable(F.T * F)
    
    Ic = ufl.variable(ufl.tr(C))
    J  = ufl.variable(ufl.det(F))

    psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J))**2
    P = ufl.diff(psi, F)

    #P = 2.0 * mu * ufl.sym(ufl.grad(u)) + lmbda * ufl.tr(ufl.sym(ufl.grad(u))) * I
    vM = ufl.sqrt(2./3. * ufl.inner(ufl.dev(P), ufl.dev(P)))
    V_vonMises = dlfn.fem.FunctionSpace(mesh, ("DG", 0)) # TODO why DG0, not CG1?
    vM_expr = dlfn.fem.Expression(vM, V_vonMises.element.interpolation_points())
    vM_Fun = dlfn.fem.Function(V_vonMises, name="vonMises")
    vM_Fun.interpolate(vM_expr)
    return vM_Fun

load_steps = 200
uval = 4 * 1e-3 / (load_steps) # m
#u_rval = 1.8 * 1e-3 / (load_steps - 1)

for n in range(1, load_steps + 2):
    print(f"------------------------------ Load step {n}: ------------------------------")
    obstacle_fun.interpolate(obstacle.eval)
    print(f"Current z_displacement u_z_presc = {obstacle.z_displacement * 1e3:.2} mm")
    obstacle_funx, obstacle_funy, obstacle_funz = obstacle_fun.split()
    num_its, converged = solver.solve(u)
    u.name = "Deflection"
    xdmf.write_function(u, n)
    u.x.scatter_forward() # TODO this has something to do with MPI stuff which I do not understand 
    assert(converged)
    print(f"Load step {n}, Iterations {num_its}, F = {get_force(u):.2} N, u_r = {obstacle.u_r * 1e3:.2} mm, max(|u_z|) = {np.max(np.abs(u.x.array.reshape(-1 , len(u))[:,2])) * 1e3:.2f} mm\n")
    #traction_foo = get_pressure(u)
    #print(f"Contact pressure = {traction_foo}")
    #t_coulomb.value = - np.abs(traction_foo) / 2

    # update displacement
    #u_D0.value[2] = n * uval
    obstacle.z_displacement = n * uval
    #obstacle.u_r = (n - 1) * u_rval
    # normal_traction.value = 0.34 * n * 8000000
    # tangential_traction.value = 0.93 * n * 8000000

    vM_Fun = get_vonMises(u)
    xdmf.write_function(vM_Fun, n)