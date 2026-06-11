# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 11:11:04 2026

@author: USER-01
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 02:30:16 2026

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


steady_state_by_intensity =[]

for j in range (0, 6, 1):
    
    laser_frequency_vals = []

    steady_state_vals = []

    for i in range(-200, 200, 5):
        # Adding the lasers
        laser_1 = library.Laser(
            frequency = transition_12.frequency() + i * 1e7,
            intensity = 10**j,
            linewidth = 4e6,
        )
        
        laser_frequency_vals.append(laser_1.frequency)
        
        laser_2 = library.Laser(
            frequency = c / 405e-9,
            intensity = 1e11,
            linewidth = 4e6,
        )
        
        lasers = [laser_1, laser_2]
        
        Omega1 = np.dot(transition_12.transitionDipoleMoment(), laser_1.electric_field_amplitude()) / hbar
        Delta1 = (laser_1.frequency - transition_12.frequency()) * 2 * np.pi
        
        rabi_period = np.pi * 2 / Omega1
        
        sigma_ion = 6e-21
        Gamma_ion = sigma_ion * laser_2.intensity / (h * laser_2.frequency)
        
        g = basis(4,0)
        e = basis(4,1)
        continuum = basis(4,2)
        
        # Hamiltonian in rotating frame
        H = (
            Omega1/2 * (g*e.dag() + e*g.dag())
            - Delta1 * e*e.dag()
            - Delta1 * continuum*continuum.dag()
        )
        
        # The jump operators
        c_ops = [
            np.sqrt(transition_12.A) * g*e.dag(),
            # np.sqrt(Gamma_ion) * c*e.dag()
        ]
        
        rho0 = g * g.dag()
        tlist = np.linspace(0, 1e-7, num=1000)
        
        result = mesolve(
            H,
            rho0,
            tlist,
            c_ops = c_ops,
            e_ops = [
                g*g.dag(),
                e*e.dag(),
                continuum*continuum.dag()
            ]
        )
        
        # taking the average of the last 1/10 of the data(rho22)
        data = result.expect[1]
        steady_state_population = np.mean(data[int(0.9*len(data)):])
        
        print(f'pho22 = {steady_state_population}')
        steady_state_vals.append(steady_state_population)
    
    steady_state_by_intensity.append(steady_state_vals)
    
#%%
lambda0 = transition_12.frequency()
# Convert wavelength scan to frequency detuning
detuning_MHz = (np.array(laser_frequency_vals) - lambda0) / 1e6

steady_state_by_intensity = np.array(steady_state_by_intensity)

fig, ax = plt.subplots(figsize=(8,5))

ax.set_xlabel('553 nm laser detuning (MHz)', fontsize=12)

ax.set_ylabel(r'Steady-state $\rho_{22}$', fontsize=12)

ax.set_title('Scanning 553 nm laser (Ba-138)', fontsize=13)

ax.plot(detuning_MHz, steady_state_by_intensity[0],
        label=r'$10^{0}\,\mathrm{W\,m^{-2}}$')
ax.plot(detuning_MHz, steady_state_by_intensity[1],
        label=r'$10^{1}\,\mathrm{W\,m^{-2}}$')
ax.plot(detuning_MHz, steady_state_by_intensity[2],
        label=r'$10^{2}\,\mathrm{W\,m^{-2}}$')
ax.plot(detuning_MHz, steady_state_by_intensity[3],
        label=r'$10^{3}\,\mathrm{W\,m^{-2}}$')
ax.plot(detuning_MHz, steady_state_by_intensity[4],
        label=r'$10^{4}\,\mathrm{W\,m^{-2}}$')
ax.plot(detuning_MHz, steady_state_by_intensity[5],
        label=r'$10^{5}\,\mathrm{W\,m^{-2}}$')

ax.grid(alpha=0.3)

ax.legend()

plt.show()

