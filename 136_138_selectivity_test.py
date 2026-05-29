# -*- coding: utf-8 -*-
"""
Created on Thu May 28 17:55:06 2026

@author: USER-01
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May 26 17:20:19 2026

@author: USER-01
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May 21 11:49:00 2026

@author: USER-01
"""

import numpy as np
import library
import matplotlib.pyplot as plt
from scipy.constants import m_p, c

# This code snippet is for plotting the ionization rate peaks of Ba136 and Ba138
# A pulsed laser of duration 1ns and intensity 1e11 is used

# Scan the laser for Ba 138
# A list to store excited population at 1ns
population_vals_138 = []

# And the corresponding laser frequencies
laser_frequency_vals_138 = []

for i in range(-100, 100, 1):
    print(i)
    # Adding the states
    state_1 = library.State(
        configuration = "6s2",
        S = 0,
        L = 0,
        J = 0,
        energy = 0
    )
    # print(ground_state.spectroscopic_name())
    
    state_2 = library.State(
        configuration = "6s6p",
        S = 0,
        L = 1,
        J = 1,
        energy = 18060.261
    )
    
    state_3 = library.State(
        configuration = "6s6d",
        S = 0,
        L = 2,
        J = 2,
        energy = 30236.826
    )
    
    states = [state_1, state_2, state_3]
    
    
    # Add the transitions
    transition_12 = library.Transition(state_1, state_2,
                                      linewidth = 20e6, dipole=None, A = 1.19e8)
    
    transition_23 = library.Transition(state_2, state_3,
                                      linewidth = 20e6, dipole=None, A = 5.00e6)
    
    transition_13 = library.Transition(state_1, state_3,
                                      linewidth = 20e6, dipole=None, A = 1e5)
    
    transitions = [transition_12, transition_23, transition_13]
    
    # Adding the lasers
    
    # Scanning the first laser, step = 10MHz
    laser_1 = library.Laser(
        frequency = transition_12.frequency() + i * 1e7,
        intensity = 0.5e3,
        linewidth = 4e6,
    )
    
    laser_frequency_vals_138.append(laser_1.frequency)
    
    laser_2 = library.Laser(
        frequency = transition_23.frequency(),
        intensity = 1e5,
        linewidth = 4e6,
    )
    
    laser_3 = library.Laser(
        frequency = c / 405e-9, 
        intensity = 1e11,
        linewidth = 4e6,
    )
    
    lasers = [laser_1, laser_2, laser_3]
    
    # Testing the solver
    barium_solver = library.RateEquationSolver(states, transitions, lasers,
                                               temperature = 500, mass = 138 * m_p)
    
    N0 = [1, 0, 0, 0]   # Initially all in ground state
    population = barium_solver.population_at_t(N0, 1e-9)
    population_vals_138.append(population)
    
    
    
# Scan the laser for Ba 136
# A list to store excited population at 1ns
population_vals_136 = []

# And the corresponding laser frequencies
laser_frequency_vals_136 = []

for i in range(-100, 100, 1):
    print(i)
    # Adding the states
    state_1 = library.State(
        configuration = "6s2",
        S = 0,
        L = 0,
        J = 0,
        energy = 0
    )
    # print(ground_state.spectroscopic_name())
    
    state_2 = library.State(
        configuration = "6s6p",
        S = 0,
        L = 1,
        J = 1,
        energy = 18060.261 + (1 / (c / 0.12802 * 1e9)) * 1e2
    )
    
    state_3 = library.State(
        configuration = "6s6d",
        S = 0,
        L = 2,
        J = 2,
        energy = 30236.826
    )
    
    states = [state_1, state_2, state_3]
    
    
    # Add the transitions
    transition_12 = library.Transition(state_1, state_2,
                                      linewidth = 20e6, dipole=None, A = 1.19e8)
    
    transition_23 = library.Transition(state_2, state_3,
                                      linewidth = 20e6, dipole=None, A = 5.00e6)
    
    transition_13 = library.Transition(state_1, state_3,
                                      linewidth = 20e6, dipole=None, A = 1e5)
    
    transitions = [transition_12, transition_23, transition_13]
    
    # Adding the lasers
    
    # Scanning the first laser, step = 10MHz
    laser_1 = library.Laser(
        frequency = transition_12.frequency() + i * 1e7,
        intensity = 0.5e3,
        linewidth = 4e6,
    )
    
    laser_frequency_vals_136.append(laser_1.frequency)
    
    laser_2 = library.Laser(
        frequency = transition_23.frequency(),
        intensity = 1e5,
        linewidth = 4e6,
    )
    
    laser_3 = library.Laser(
        frequency = c / 405e-9, 
        intensity = 1e11,
        linewidth = 4e6,
    )
    
    lasers = [laser_1, laser_2, laser_3]
    
    # Testing the solver
    barium_solver = library.RateEquationSolver(states, transitions, lasers,
                                               temperature = 500, mass = 136 * m_p)
    
    N0 = [1, 0, 0, 0]   # Initially all in ground state
    population = barium_solver.population_at_t(N0, 1e-9)
    population_vals_136.append(population)
    

fig, ax = plt.subplots(figsize=(8, 5))

lambda0_138 = transition_12.frequency()

lambda0_136 = transition_12.frequency() + 0.12802 * 1e9

# Convert wavelength scan to frequency detuning
detuning_MHz_138 = (np.array(laser_frequency_vals_138) - lambda0_138) / 1e6

population_vals_138 = np.array(population_vals_138)

# Same for Ba 136
detuning_MHz_136 = (np.array(laser_frequency_vals_136) - lambda0_136) / 1e6

population_vals_136 = np.array(population_vals_136)

fig, ax = plt.subplots(figsize=(8,5))



ax.plot(detuning_MHz_138, population_vals_138[:,3],
        label='Ba138')

ax.plot(detuning_MHz_136, population_vals_136[:,3],
        label='Ba136')

ax.set_xlabel('553 nm laser detuning (MHz)', fontsize=12)

ax.set_ylabel('Population fraction at 1 ns', fontsize=12)

ax.set_title('Ba-136 and Ba-138 peak test', fontsize=13)

ax.grid(alpha=0.3)

ax.legend()

plt.show()

'''
population_vals = np.array(population_vals)
ax.plot(laser_frequency_vals, population_vals[:,0], label = '$N_1$ ground (6s²)')
ax.plot(laser_frequency_vals, population_vals[:,1], label = '$N_2$ (6s6p ¹P₁)')
ax.plot(laser_frequency_vals, population_vals[:,2], label = '$N_3$ (6s6d ¹D₂)')
ax.plot(laser_frequency_vals, population_vals[:,3], label = '$N_4$ ion')
ax.set_xlabel('553nm laser wavelength', fontsize=12)
ax.set_ylabel('Ion fraction at 1ns', fontsize=12)
ax.set_title('Scanning the 553nm laser_test_Ba138', fontsize=13)
plt.legend()
plt.show()
'''
    
