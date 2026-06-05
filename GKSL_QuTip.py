# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 16:45:27 2026

@author: USER-01
"""

import numpy as np
import library
from qutip import *
from scipy.constants import m_p, c, hbar, epsilon_0
import matplotlib.pyplot as plt

print(qutip.__version__)

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
                                  linewidth = 20e6, lifetime = None, A = 1.19e8)

transition_23 = library.Transition(state_2, state_3,
                                  linewidth = 20e6, lifetime=None, A = 5.00e6)

transition_13 = library.Transition(state_1, state_3,
                                  linewidth = 20e6, lifetime=None, A = 1e5)

transitions = [transition_12, transition_23, transition_13]

# Adding the lasers
# Testing the on resonance condition
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


Omega1 = np.dot(transition_12.transitionDipoleMoment(), laser_1.electric_field_amplitude()) / hbar
Omega2 = np.dot(transition_23.transitionDipoleMoment(), laser_2.electric_field_amplitude()) / hbar

Delta1 = (laser_1.frequency - transition_12.frequency()) * 2 * np.pi
Delta2 = (laser_2.frequency - transition_23.frequency()) * 2 * np.pi

g = basis(4,0)
e = basis(4,1)
r = basis(4,2)
c = basis(4,3)

H = (
    Omega1/2 * (g*e.dag() + e*g.dag())
    + Omega2/2 * (e*r.dag() + r*e.dag())
    - Delta1 * e*e.dag()
    - (Delta1+Delta2) * r*r.dag()
)

# The linewidth of the final transition is a guess!
c_ops = [
    np.sqrt(transition_12.A) * g*e.dag(),
    np.sqrt(transition_23.A) * e*r.dag(),
    np.sqrt(1e8) * c*r.dag()
]

rho0 = g * g.dag()
tlist = np.linspace(0, 1e-9, num=200)

result = mesolve(
    H,
    rho0,
    tlist,
    c_ops = c_ops,
    e_ops = [
        g*g.dag(),
        e*e.dag(),
        r*r.dag(),
        c*c.dag()
    ]
)

fig, ax = plt.subplots()

ax.plot(tlist, result.expect[0]);

ax.plot(tlist, result.expect[1]);

ax.plot(tlist, result.expect[2]);

ax.plot(tlist, result.expect[3]);

ax.set_xlabel('Time');

ax.set_ylabel('Population fraction');

ax.legend([r'$\rho_{11}$', r'$\rho_{22}$', r'$\rho_{33}$', 'Ion population'])

plt.show(fig)

rho_ss = steadystate(H, c_ops)

print(rho_ss)