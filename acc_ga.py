# Python Accelerometer GA
# @ruiesteves

# Imports
import acc_geo
import acc_disp
import acc_elec
import acc_modal
import acc_damping
import numpy as np
import math as math
import random as rand
import copy
import time

# Initial parameters of the device to be optimized
scale = 1e-6
suspension_beam_width = 350*scale
proof_mass_length = 2400*scale

initial = [suspension_beam_width,proof_mass_length]

# Classes
class GA_device:

    def __init__(self,id):
        self.list_parameters = []
        self.id = id

    def calc_sensitivity(self):
        list = self.list_parameters
        self.sensitivity = acc_elec.main(list[0],list[1])

    def calc_freq(self):
        list = self.list_parameters
        self.freq = acc_modal.main(list[0],list[1])

    def calc_qfactor(self):
        list = self.list_parameters
        self.qfactor = acc_damping.q_factor(list[0],list[1],self.freq)

    def calc_fom(self):
        list = self.list_parameters
        try:
            acc_geo.build(list[0],list[1])
            self.calc_sensitivity()
            self.calc_freq()
            self.calc_qfactor()
            self.fom = self.freq * self.sensitivity * self.qfactor * (1/1000)
        except:
            print("Geometry became invalid for device",self.id)
            self.fom = 0


class GA:   # GA class, initiated with a list of devices, a list of mutation chances and a list of mutation relative size

    def __init__(self,list_devices,mutation_chance,mutation_size):
        self.list_devices = list_devices            # Must be a list of GA_devices
        self.mutation_chance = mutation_chance      # A list, with different (or not) mutation chances for each parameter
        self.mutation_size = mutation_size          # Same as above, this time for mutation_sizes (IMPORTANT to check)

    def mutate(self,dev):                           # The mutation function, mutating the parameters according to their mutation chance and size
        le = len(dev.list_parameters)
        rand.seed(time.time())
        for i in range(le):
            if rand.uniform(0,1) < self.mutation_chance[i]:
                if rand.uniform(0,1) < 0.5:
                    dev.list_parameters[i] = dev.list_parameters[i] + dev.list_parameters[i]*self.mutation_size[i]
                else:
                    dev.list_parameters[i] = dev.list_parameters[i] - dev.list_parameters[i]*self.mutation_size[i]

    def reproduce(self,top_25):
        new_population = []

        for dev in top_25:                                      # Passing the best 25 devices to the next generation
            new_population.append(dev)

        for dev in top_25:                                      # Copying and mutating the best 25 devices to the next generation
            new_dev = copy.deepcopy(dev)
            self.mutate(new_dev)
            new_population.append(new_dev)

        for i in range(len(self.list_devices)//2):              # Randomly mutating and passing half of the population to the next generation
            new_dev_r = copy.deepcopy(self.list_devices[i])
            self.mutate(new_dev_r)
            new_population.append(new_dev_r)

        return new_population


    def one_generation(self):

        for dev in self.list_devices:
            dev.calc_fom()

        le = len(self.list_devices)
        scores = [self.list_devices[i].fom for i in range(le)]
        max = np.amax(scores)
        print(scores)
        print(max)

        top_25_index = list(np.argsort(scores))[3*(le//4):le]
        top_25 = [self.list_devices[i] for i in top_25_index][::-1]

        self.list_devices = self.reproduce(top_25)





# Script
print("Genetic algorithm optimization for MEMS accelerometer")
num_pop = int(input("Size of the population: "))
num_gen = int(input("Number of generations: "))

initial_pop = []
for i in range(num_pop):
    initial_pop.append(GA_device(i))
    for par in range(len(initial)):
        initial_pop[i].list_parameters.append(initial[par])
    initial_pop[i].calc_fom()

init_ga = GA(initial_pop,[0.6,0.6],[0.13,0.0143])

for i in range(num_gen):
    init_ga.one_generation()
    for dev in init_ga.list_devices:
        print("\n","For Device number",dev.id,":")
        print(dev.sensitivity*1e3,"mV/g")
        print(dev.freq,"Hz")
        print(dev.qfactor,"Q-factor")
        print(dev.fom,"FOM")

max_score = 0


for dev in init_ga.list_devices:
    if dev.fom > max_score:
        max_score = dev.fom

print("Maximum FOM:",max_score)
