# This file shows how to use the features of the class Membrane
from __future__ import print_function
import espressomd
from Membrane import Membrane

#For Visualization
from espressomd.visualization_opengl import *
from threading import Thread
from time import sleep


import numpy as np

# System parameters
#############################################################
numLipids   = 500
skin 		= 0.4
time_step	= 0.001
eq_tstep 	= 0.0001
temperature = 1.1
box_l 		= 17.

warm_steps 	= 100
warm_n_time = 2000
min_dist 	= 0.7

# integration
sampling_interval 		= 100
equilibration_interval 	= 1000
sampling_iterations 	= 100
equilibration_iterations= 5

# Interaction parameters (Lennard Jones)
#############################################################

lj_eps = 1.0
lj_sig = 1.0
lj_cut = pow(2., 1. / 6) * lj_sig
lj_cap = 2.
lj_cut_mixed = pow(2., 1. / 6) * lj_sig

# Set Up
#############################################################
system = espressomd.System()

system.time_step 		= time_step
system.cell_system.skin = skin
system.box_l 			= [box_l, box_l, box_l]

# Membrane
#############################################################
membrane = Membrane(system,numLipids)
membrane.setOrientation("bilayer")

# Non bonded Interactions between the beads
#############################################################
system.non_bonded_inter[0, 0].lennard_jones.set_params(
        epsilon=lj_eps, sigma=0.95 * lj_sig,
        cutoff=lj_cut, shift=1. / 4)

system.non_bonded_inter[1, 1].lennard_jones.set_params(
        epsilon=lj_eps, sigma=lj_sig,
        cutoff=lj_cut, shift=1. / 4)

system.non_bonded_inter[0, 1].lennard_jones.set_params(
        epsilon=lj_eps, sigma=0.95 * lj_sig,
        cutoff=lj_cut_mixed, shift=1. / 4)

# Attractive Tail-Tail
system.non_bonded_inter[1, 1].lennard_jones_cos2.set_params(
        epsilon=lj_eps, sigma=lj_sig,
        width=1.6 * lj_sig, offset=0.)



print("\nThe position of the beads in 97th lipid: ", membrane.lipid[97].pos)
print("\nThe particle Id's for the beads in 97th lipid: ",membrane.lipid[97].partId)
print("\nThe type of particles that make the 97th lipid: ",membrane.lipid[97].type)
