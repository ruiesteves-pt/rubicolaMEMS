# MEMS accelerometer squeeze film damping calculation script
# @ruiesteves

import math

epsilon0 = 8.85e-12
scale = 1e-6
mu = 1.86e-5            # the mean viscosity of the medium
lamb = 0.067e-6         # mean free path
small_gap = 22*scale
thickness = 320*scale
rho = 2329

def q_factor(suspension_beam_width,proof_mass_length,sense_frequency):
    A = proof_mass_length**2
    mass_sense = A * thickness * rho
    Kn = lamb/small_gap
    mu_eff = mu/(1+9.638*Kn**1.159)
    Pa = 101.3e3
    c = 1
    squeeze_number = (12*mu_eff*2*math.pi*L_sense**2)/(Pa*small_gap**2)
    sum = 0
    for m in range(1,10,2):
        for n in (1,10,2):
            sum = sum + (m**2 + c**2 * n**2)/((m*n)**2 * ((m**2 + c**2 * n**2)**2 + (squeeze_number**2 / math.pi**4)))

    F_damping = ((64*squeeze_number*Pa*A)/(math.pi**6 * small_gap)) * sum
    c_sense = F_damping
    q_factor = mass_sense * sense_frequency*2*math.pi / c_sense

    return q_factor
