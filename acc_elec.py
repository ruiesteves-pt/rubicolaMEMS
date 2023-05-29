# MEMS accelerometer electrical domain simulation
# @ruiesteves

# Imports
import acc_disp



# Constants
e0 = 8.85e-12           # Permitivity of free space
er = 1                  # Relative permitivity of dielectric, in this case air
small_gap = 22*scale    # Distance between electrodes and proof mass


# Functions

def main(suspension_beam_width,proof_mass_length):
    disp = acc_disp.disp(suspension_beam_width,proof_mass_length)
    A = (proof_mass_length)**2
    def top_capacitance():
        C = (e0*er*A)/(small_gap + disp)
        return C

    def bot_capacitance():
        C = (e0*er*A)/(small_gap - disp)
        return C

    def capacitance_total():
        c_total = bot_capacitance() - top_capacitance()
        print(c_total*1e15,"fF")
        return c_total

    def c2v():
        v = 2 * capacitance_total() * 2.5 * (1/300e-15)
        print("Output voltage:",v*1e3,"mV")
        return v

    voltage = c2v()
    return voltage
