# -*- coding: utf-8 -*-
"""
Created on Fri May 29 10:46:58 2026

@author: USER-01
"""

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

# This code snippet is for estimating isotope selectivity at different temperatures
# A pulsed laser of duration 1ns and intensity 1e11 is used

temperature_vals = []

# Lists to store isotope fractions(divided by Ba138 population)
Ba_136_fraction = []

for i in range (0, 100, 1):
    temperature = i * 10
    temperature_vals.append(temperature)

    # Finding out Ba138 population at 1ns
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
    laser_1_138 = library.Laser(frequency = transition_12_138.frequency(), intensity = 0.5e3, linewidth = 4e6)
    laser_2_138 = library.Laser(frequency = transition_23_138.frequency(), intensity = 1e5,linewidth = 4e6)
    laser_3_138 = library.Laser(frequency = c / 405e-9, intensity = 1e11, linewidth = 4e6)
    lasers_138 = [laser_1_138, laser_2_138, laser_3_138]
        
    barium_solver_138 = library.RateEquationSolver(states_138, transitions_138, lasers_138,
                                               temperature = temperature, mass = 138 * m_p)
    
    N0_138 = [1, 0, 0, 0]   # Initially all in ground state
    population_138 = barium_solver_138.population_at_t(N0_138, 1e-9)

    # Finding out the population of Ba136 at Ba138 resonance
    population_vals_136 = []
    state_1_136 = library.State(configuration = "6s2", S = 0, L = 0, J = 0, energy = 0)
    state_2_136 = library.State(configuration = "6s6p", S = 0, L = 1, J = 1, energy = 18060.261 + (1 / (c / 0.12802 * 1e9)) * 1e2)
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
    lasers_136 = lasers_138
    
    barium_solver_136 = library.RateEquationSolver(states_136, transitions_136, lasers_136,
                                               temperature = temperature, mass = 136 * m_p)
    
    N0_136 = [1, 0, 0, 0]   # Initially all in ground state
    population_136 = barium_solver_136.population_at_t(N0_136, 1e-9)
    
    
    Ba_136_fraction.append(population_136[3] / population_138[3])



fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(temperature_vals, Ba_136_fraction,
        label='Ba136')

ax.set_xlabel('Temperature(K)', fontsize=12)
ax.set_ylabel('Ba isotope fractions', fontsize=12)

ax.set_title('Temperature dependence of isotope selectivity', fontsize=13)
ax.grid(alpha=0.3)
ax.legend()

plt.show()


    
