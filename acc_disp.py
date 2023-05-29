# MEMS accelerometer displacement simulation script
# @ruiesteves

# Imports
from __future__ import print_function
from dolfin import *

# Material constants
E = Constant(170e9)
nu = Constant(0.28)
rho = 2329
mu = E/2/(1+nu)
lmbda = E*nu/(1+nu)/(1-2*nu)

# Mesh
mesh = Mesh('accelerometer.xml')


def disp(suspension_beam_width,proof_mass_length):

    # Constants
    cl = 175
    proof_mass_cl = 150
    scale = 1e-6
    suspension_beam_length = 3300*scale
    beam_thickness = 69*scale
    small_beam_length = 500*scale
    proof_mass_thickness = 320*scale
    beam_dist = 500*scale
    beam_l = 122.5 * scale
    beam_h = 177.5 * scale
    beam_dist2 = beam_dist + suspension_beam_length
    beam_dist3 = proof_mass_length + suspension_beam_width + small_beam_length
    beam_dist_final = beam_dist2 - beam_dist3
    beam_lower = (proof_mass_thickness - beam_thickness)/2
    beam_to_mass = (beam_dist_final+suspension_beam_width+small_beam_length + proof_mass_length) - suspension_beam_length
    anchor_top = beam_dist_final+suspension_beam_width+small_beam_length+beam_to_mass+suspension_beam_length
    volume_PM = proof_mass_length**2 * proof_mass_thickness

    # Strain operator
    def eps(v):
        return sym(grad(v))

    # Stress tensor
    def sigma(v):
        return lmbda*tr(eps(v))*Identity(3) + 2.0*mu*eps(v)


    # Boundary
    def left(x,on_boundary):
        return near(x[0],0.)

    def bottom(x,on_boundary):
        return near(x[1],0.)

    def top(x,on_boundary):
        return near(x[1],anchor_top)

    def right(x,on_boundary):
        return near(x[0],anchor_top)

    def main():
        rho_g = 9.8*(rho)
        print("Volume:",volume_PM)
        print("Force:",rho_g)
        f = Constant((0.,0.,rho_g))
        T = Constant((0,0,0))
        V = VectorFunctionSpace(mesh,'Lagrange',degree=3)
        du = TrialFunction(V)
        u_ = TestFunction(V)
        a = inner(sigma(du),eps(u_))*dx
        l = dot(f,u_)*dx

        bc = [DirichletBC(V, Constant((0.,0.,0.)),left),
        DirichletBC(V, Constant((0.,0.,0.)),right),
        DirichletBC(V, Constant((0.,0.,0.)),top),
        DirichletBC(V, Constant((0.,0.,0.)),bottom)]
        u = Function(V, name='Displacement')
        solve(a == l, u, bc)
        z_disp = u(beam_dist_final+suspension_beam_width+small_beam_length+proof_mass_length/2,beam_dist_final+suspension_beam_width+small_beam_length+proof_mass_length/2,proof_mass_thickness/2)

        # Set up file for exporting results
        file_results = XDMFFile("acc_displacement.xdmf")
        file_results.parameters["flush_output"] = True
        file_results.parameters["functions_share_mesh"] = True
        file_results.write(u,0)

        print(z_disp[2]*1e6,"um")
        return z_disp[2]

    return main()
