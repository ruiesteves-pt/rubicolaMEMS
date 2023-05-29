# MEMS accelerometer modal analysis script
# @ruiesteves

# Imports
from fenics import *
import numpy as np

# Material constants
E = Constant(170e9)
nu = Constant(0.28)
rho = 2329
mu = E/2/(1+nu)
lmbda = E*nu/(1+nu)/(1-2*nu)

# Meshing
mesh = Mesh('accelerometer.xml')


def main(suspension_beam_width,proof_mass_length):

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

    # Functions
    def eps(v):
        return sym(grad(v))

    def sigma(v):
        return lmbda*tr(eps(v))*Identity(3) + 2.0*mu*eps(v)

    # Function Space
    V = VectorFunctionSpace(mesh,'Lagrange',degree=3)
    u_ = TrialFunction(V)
    du = TestFunction(V)

    # Boundary
    def left(x,on_boundary):
        return near(x[0],0.)

    def bottom(x,on_boundary):
        return near(x[1],0.)

    def top(x,on_boundary):
        return near(x[1],anchor_top)

    def right(x,on_boundary):
        return near(x[0],anchor_top)

    bc = [DirichletBC(V, Constant((0.,0.,0.)),left),
    DirichletBC(V, Constant((0.,0.,0.)),right),
    DirichletBC(V, Constant((0.,0.,0.)),top),
    DirichletBC(V, Constant((0.,0.,0.)),bottom)]

    # Matrices
    k_form = inner(sigma(du),eps(u_))*dx
    l_form = Constant(1.)*u_[0]*dx
    K = PETScMatrix()
    b = PETScVector()
    assemble_system(k_form,l_form,bc,A_tensor=K,b_tensor=b)

    m_form = rho*dot(du,u_)*dx
    M = PETScMatrix()
    assemble(m_form, tensor=M)

    # Eigenvalues/Eigensolver
    eigensolver = SLEPcEigenSolver(K,M)
    eigensolver.parameters['problem_type'] = 'gen_hermitian'
    #eigensolver.parameters['spectrum'] = 'smallest real'
    eigensolver.parameters['spectral_transform'] = 'shift-and-invert'
    eigensolver.parameters['spectral_shift'] = 0.
    N_eig = 2
    eigensolver.solve(N_eig)
    #print (eigensolver.parameters.str(True))

    # Export results
    file_results = XDMFFile('acc_modal_analysis.xdmf')
    file_results.parameters['flush_output'] = True
    file_results.parameters['functions_share_mesh'] = True

    r1,c1,rx1,cx1 = eigensolver.get_eigenpair(0)
    u = Function(V)
    u.vector()[:] = rx1
    file_results.write(u,0)

    # Extraction
    for i in range(N_eig):
        r,c,rx,cx = eigensolver.get_eigenpair(i)
        freq = sqrt(r)/2/pi
        print('Mode:',i,'   ','Freq:',freq,'[Hz]')


    freq_final = sqrt(r1)/2/pi
    return freq_final
