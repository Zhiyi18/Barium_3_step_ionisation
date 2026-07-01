# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 13:56:15 2026

@author: USER-01
"""

import numpy as np
import library
from qutip import *
from scipy.constants import m_p, c, h, hbar, epsilon_0
import matplotlib.pyplot as plt

# Adding the states
state_1 = library.State(
    configuration = "6s2",
    S = 0,
    L = 0,
    J = 0,
    energy = 0
)

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

J_states = [state_1, state_2, state_3]

# Add the hyperfine states
hyperfine_states = []
for states in J_states:
    hyperfine_state = states.hyperfine_states(I=3/2, A=-109.5, B=50.09)
    hyperfine_states.append(hyperfine_state)
    
# print(hyperfine_states[2][3].energy)

'''
Now the hyperfine states are stored like [[hyperfine states of the first J state], [the second J state], ...]
Then we add the transitions without hyperfine splitting to find the 'base' laser frequency
'''   
transition_12 = library.Transition(state_1, state_2,
                                  linewidth = 20e6, lifetime = None, A = 1.19e8)

transition_23 = library.Transition(state_2, state_3,
                                  linewidth = 20e6, lifetime=None, A = 5.00e6)

transition_13 = library.Transition(state_1, state_3,
                                  linewidth = 20e6, lifetime=None, A = 1e5)


laser_1 = library.Laser(
    frequency = transition_12.frequency(),
    intensity = 1e5,
    linewidth = 4e6,
)

laser_2 = library.Laser(
    frequency = transition_23.frequency(),
    intensity = 1e5,
    linewidth = 4e6,
)

laser_3 = library.Laser(
    frequency = c / 405e-9,
    intensity = 1e10,
    linewidth = 4e6,
)


'''
What to do then:
    calculate the transition diople moments using the clebsch-gordan coefficient function from qutip
    add the lasers
    calculate the rabi frequencies
    define the detunings like Delta1 = (laser_1.frequency - transition_12.frequency()) * 2 * np.pi
    construct the 8*8 hamiltonian(for 8 hyperfine states)
    construct the jump operators
    pass into qutip GKSL solver and plot the populations in 8 levels
'''





