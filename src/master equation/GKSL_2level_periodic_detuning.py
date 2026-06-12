# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 11:30:07 2026

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
# print(ground_state.spectroscopic_name())

state_2 = library.State(
    configuration = "6s6p",
    S = 0,
    L = 1,
    J = 1,
    energy = 18060.261
)

states = [state_1, state_2]


# Add the transitions
transition_12 = library.Transition(state_1, state_2,
                                  linewidth = 20e6, lifetime=None, A = 1.19e8)

transitions = [transition_12]

# Adding the lasers
laser_1 = library.Laser(
    frequency = transition_12.frequency(),
    intensity = 1e5,
    linewidth = 4e6,
)

laser_2 = library.Laser(
    frequency = c / 405e-9,
    intensity = 1e11,
    linewidth = 4e6,
)

lasers = [laser_1, laser_2]

Omega1 = np.dot(transition_12.transitionDipoleMoment(), laser_1.electric_field_amplitude()) / hbar
Delta1 = (laser_1.frequency - transition_12.frequency()) * 2 * np.pi

rabi_period = np.pi * 2 / Omega1
print(f'period = {rabi_period}')

sigma_ion = 6e-21
Gamma_ion = sigma_ion * laser_2.intensity / (h * laser_2.frequency)

g = basis(4,0)
e = basis(4,1)
c = basis(4,2)

# Hamiltonian in rotating frame
H0 = (Omega1/2) * (g*e.dag() + e*g.dag())

H_detuning = -(e*e.dag() + c*c.dag())

# Modulating the detuning
Delta0 = 0
# A = 20e8 * 2*np.pi      # modulation amplitude (rad/s)
# omega_m = 2*np.pi*5e6   # modulation frequency
A = Omega1
omega_m = Omega1


# QuTip takes in H0 + H_detuning * Delta_t in this way
def Delta_t(t, args):
    return Delta0 + A * np.cos(omega_m*t)

H = [
    H0,
    [H_detuning, Delta_t]
]

# The jump operators
c_ops = [
    np.sqrt(transition_12.A) * g*e.dag(),
    # np.sqrt(Gamma_ion) * c*e.dag()
]

rho0 = g * g.dag()
tlist = np.linspace(0, 1e-9, num=1000)

result = mesolve(
    H,
    rho0,
    tlist,
    c_ops = c_ops,
    e_ops = [
        g*g.dag(),
        e*e.dag(),
        c*c.dag()
    ]
)


total_population = result.expect[0] + result.expect[1] + result.expect[2]

fig, ax = plt.subplots()

ax.plot(tlist, result.expect[0], lw = 0.5);

ax.plot(tlist, result.expect[1], lw = 0.5);

ax.plot(tlist, result.expect[2], lw = 1);

ax.plot(tlist, total_population, lw = 0.3);

ax.set_xlabel('Time');

ax.set_ylabel('Population fraction');

ax.set_title('Ba 2-step ionisation — master equation model', fontsize=13)

ax.legend([r'$\rho_{11}$ (6s²)', 
           r'$\rho_{22}$ (6s6p ¹P₁)', 
           'Ion population',
           'Total population'])
