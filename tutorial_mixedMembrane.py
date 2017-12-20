# This file shows how to use the features of the class Membrane
from __future__ import print_function
import espressomd
from Lipid import Lipid
from Membrane import Membrane

#For Visualization
from espressomd.visualization_opengl import *
from threading import Thread
from time import sleep

import numpy as np


# System parameters
#############################################################
numLipids 	= 320
skin 		= 0.4
time_step 	= 0.001
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
visualizer = openGLLive(system)

system.thermostat.set_langevin(kT=temperature, gamma=1.0)


def main():
	system.time_step = time_step
	system.cell_system.skin = skin
	system.box_l = [box_l, box_l, box_l]
	lj_cap = 2.5
	
	#Membrane
	############################################################
	lipid1   = Lipid(system, sigma = 0.95,lipidType = {"Head": 0, "Mid": 1, "Tail": 1})
	lipid2   = Lipid(system, sigma = 1.05,lipidType = {"Head": 0, "Mid": 2, "Tail": 1})
	membrane = Membrane(system, numLipids)
	membrane.setOrientation('mixedbilayer',lipid1=lipid1,lipid2=lipid2)

	# Non bonded Interactions between the beads (Have to turn off the Mid-Mid cosine interaction in two types of lipids)
	#############################################################
	system.non_bonded_inter[0, 0].lennard_jones.set_params(
			epsilon=lj_eps, sigma=0.95 * lj_sig,
			cutoff=lj_cut, shift=1. / 4)

	system.non_bonded_inter[0, 1].lennard_jones.set_params(
			epsilon=lj_eps, sigma=0.95 * lj_sig,
			cutoff=lj_cut_mixed, shift=1. / 4)

	system.non_bonded_inter[0, 2].lennard_jones.set_params(
			epsilon=lj_eps, sigma=0.95 * lj_sig,
			cutoff=lj_cut_mixed, shift=1. / 4)

	system.non_bonded_inter[1, 1].lennard_jones.set_params(
			epsilon=lj_eps, sigma=lj_sig,
			cutoff=lj_cut, shift=1. / 4)

	system.non_bonded_inter[2, 1].lennard_jones.set_params(
			epsilon=lj_eps, sigma=lj_sig,
			cutoff=lj_cut, shift=1. / 4)

	system.non_bonded_inter[2, 2].lennard_jones.set_params(
			epsilon=lj_eps, sigma=lj_sig,
			cutoff=lj_cut, shift=1. / 4)

	# Attractive Tail-Tail
	system.non_bonded_inter[1, 1].lennard_jones_cos2.set_params(
			epsilon=lj_eps, sigma=lj_sig,
			width=1.6 * lj_sig, offset=0.)

	system.non_bonded_inter[2, 2].lennard_jones_cos2.set_params(
			epsilon=lj_eps, sigma=lj_sig,
			width=1.6 * lj_sig, offset=0.)

	# Important for warmup
	system.force_cap = lj_cap

	#############################################################
	#  Warmup Integration                                       #
	#############################################################

	print("""
	Start warmup integration:
	At maximum {} times {} steps
	Stop if minimal distance is larger than {}
	""".strip().format(warm_n_time, warm_steps, min_dist))

	i = 0
	act_min_dist = system.analysis.min_dist()
	while i < warm_n_time and act_min_dist < min_dist:
		system.integrator.run(warm_steps)
		act_min_dist = system.analysis.min_dist()
		print("run {} at time = {} (LJ cap= {} ) min dist = {}".strip().format(
			i, system.time, lj_cap, act_min_dist))
		i += 1
		lj_cap += 1.0
		system.force_cap = lj_cap

	system.force_cap = 0.

	print("\nWarm up finished\n")

	while(True):
		system.integrator.run(100)
		visualizer.update()


#Start simulation in seperate thread
t = Thread(target=main)
t.daemon = True
t.start()

visualizer.start()

