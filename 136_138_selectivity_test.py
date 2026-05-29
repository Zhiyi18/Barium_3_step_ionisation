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

# A list to store excited population at 1ns
population_vals_138 = []

# And the corresponding laser frequencies
laser_frequency_vals_138 = []

MHz_to_cm = 1e6 / (c * 100)
shift_136 = -128.9 * MHz_to_cm

temperature = 500


# Scan the 553nm laser for Ba138, step = 10MHz
for i in range(-100, 100, 1):
    # Add the states
    state_1_138 = library.State(configuration = "6s2", S = 0, L = 0, J = 0, energy = 0)
    state_2_138 = library.State(configuration = "6s6p", S = 0, L = 1, J = 1, energy = 18060.261)  
    state_3_138 = library.State(configuration = "6s6d", S = 0, L = 2, J = 2, energy = 30236.82)
    states_138 = [state_1_138, state_2_138, state_3_138]
    
    
    # Add the transitions
    transition_12_138 = library.Transition(state_1_138, state_2_138,
                                      linewidth = 20e6, dipole=None, A = 1.19e8)
    transition_23_138 = library.Transition(state_2_138, state_3_138,
                                      linewidth = 20e6, dipole=None, A = 5.00e6)
    transition_13_138 = library.Transition(state_1_138, state_3_138,
                                      linewidth = 20e6, dipole=None, A = 1e5)
    transitions_138 = [transition_12_138, transition_23_138, transition_13_138]
    
    # Add the lasers
    laser_1_138 = library.Laser(frequency = transition_12_138.frequency() + i * 1e7, intensity = 0.5e3, linewidth = 4e6)
    laser_2_138 = library.Laser(frequency = transition_23_138.frequency(), intensity = 1e5,linewidth = 4e6)
    laser_3_138 = library.Laser(frequency = c / 405e-9, intensity = 1e11, linewidth = 4e6)
    lasers_138 = [laser_1_138, laser_2_138, laser_3_138]
    
    laser_frequency_vals_138.append(laser_1_138.frequency)
        
    barium_solver_138 = library.RateEquationSolver(states_138, transitions_138, lasers_138,
                                               temperature = temperature, mass = 138 * m_p)
    
    N0 = [1, 0, 0, 0]   # Initially all in ground state
    population_138 = barium_solver_138.population_at_t(N0, 1e-9)
    population_vals_138.append(population_138[3])
    
    
# Scan the laser for Ba 136
population_vals_136 = []

laser_frequency_vals_136 = []

for i in range(-100, 100, 1):
    state_1_136 = library.State(configuration = "6s2", S = 0, L = 0, J = 0, energy = 0)
    state_2_136 = library.State(configuration = "6s6p", S = 0, L = 1, J = 1, energy = 18060.261 + shift_136)
    state_3_136 = library.State(configuration = "6s6d", S = 0, L = 2, J = 2, energy = 30236.826)
    states_136 = [state_1_136, state_2_136, state_3_136]
    
    # Add the transitions
    transition_12_136 = library.Transition(state_1_136, state_2_136,
                                      linewidth = 20e6, dipole=None, A = 1.19e8)
    transition_23_136 = library.Transition(state_2_136, state_3_136,
                                      linewidth = 20e6, dipole=None, A = 5.00e6)
    transition_13_136 = library.Transition(state_1_136, state_3_136,
                                      linewidth = 20e6, dipole=None, A = 1e5)
    
    transitions_136 = [transition_12_136, transition_23_136, transition_13_136]
    
    # Shine Ba138 laser on Ba136
    laser_1_136 = library.Laser(frequency = transition_12_136.frequency() + i * 1e7, intensity = 0.5e3, linewidth = 4e6)
    laser_2_136 = library.Laser(frequency = transition_23_136.frequency(), intensity = 1e5,linewidth = 4e6)
    laser_3_136 = library.Laser(frequency = c / 405e-9, intensity = 1e11, linewidth = 4e6)
    lasers_136 = [laser_1_136, laser_2_136, laser_3_136]
    
    laser_frequency_vals_136.append(laser_1_136.frequency)
    
    barium_solver_136 = library.RateEquationSolver(states_136, transitions_136, lasers_136,
                                               temperature = temperature, mass = 136 * m_p)
    
    N0 = [1, 0, 0, 0]   # Initially all in ground state
    population_136 = barium_solver_136.population_at_t(N0, 1e-9)
    population_vals_136.append(population_136[3])
    

fig, ax = plt.subplots(figsize=(8, 5))

lambda0_138 = transition_12_138.frequency()

lambda0_136 = transition_12_136.frequency()

# Convert wavelength scan to frequency detuning
detuning_MHz_138 = np.array(laser_frequency_vals_138) / 1e6


population_vals_138 = np.array(population_vals_138)

# Same for Ba 136
detuning_MHz_136 = np.array(laser_frequency_vals_136) / 1e6

population_vals_136 = np.array(population_vals_136)

fig, ax = plt.subplots(figsize=(8,5))



ax.plot(detuning_MHz_138, population_vals_138,
        label='Ba138')

ax.plot(detuning_MHz_136, population_vals_136,
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
    
