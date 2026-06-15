# -*- coding: utf-8 -*-
"""
Created on Thu May 21 11:49:00 2026

@author: USER-01
"""

import numpy as np
import library
import matplotlib.pyplot as plt
from scipy.constants import m_p, c

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
                                  linewidth = 20e6, lifetime=None, A = 1.19e8)

transition_23 = library.Transition(state_2, state_3,
                                  linewidth = 20e6, lifetime=None, A = 5.00e6)

transition_13 = library.Transition(state_1, state_3,
                                  linewidth = 20e6, lifetime=None, A = 1e5)

transitions = [transition_12, transition_23, transition_13]

# Adding the lasers
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
    intensity = 1e11,
    linewidth = 4e6,
)

lasers = [laser_1, laser_2, laser_3]

# Testing the solver
barium_solver = library.RateEquationSolver(states, transitions, lasers,
                                           temperature = 500, mass = 138 * m_p)
t_vals = np.linspace(0, 1e-9, num = 1000, endpoint=True)

N0 = [1, 0, 0, 0]   # Initially all in ground state
population_vals = barium_solver.solve(N0, t_vals)

#print(population_vals)

total_population = population_vals[:,0] + population_vals[:,1] + population_vals[:,2] + population_vals[:,3]
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(t_vals, population_vals[:,0], label = '$N_1$ ground (6s²)')
ax.plot(t_vals, population_vals[:,1], label = '$N_2$ (6s6p ¹P₁)')
ax.plot(t_vals, population_vals[:,2], label = '$N_3$ (6s6d ¹D₂)')
ax.plot(t_vals, population_vals[:,3], label = '$N_4$ ion')
ax.plot(t_vals, total_population, label = 'Total population')

ax.set_xlabel('Time (s)', fontsize=12)
ax.set_ylabel('Fractional population', fontsize=12)
# ax.set_title('Transient behavior on $10^{-11}$s timescale', fontsize=13)
ax.set_title('Ba 3-step ionisation — rate equation model', fontsize=13)
# ax.set_title(r'Ionization cross-section =  $6^{-21}$', fontsize=13)
ax.legend(fontsize=11)

plt.show()

print(population_vals[:,0] + population_vals[:,1] + population_vals[:,2] + population_vals[:,3] )
