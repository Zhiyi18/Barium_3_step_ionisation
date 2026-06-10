# -*- coding: utf-8 -*-
"""
Spyderエディタ

これは一時的なスクリプトファイルです。
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, hbar, epsilon_0, c, k
from scipy.linalg import expm
from scipy.special import voigt_profile

class State:
    def __init__(self, configuration, S, L, J, energy, 
                 parity=None, mJ=None):
        
        self.configuration = configuration
        
        # Quantum numbers
        self.S = S
        self.L = L
        self.J = J
        self.mJ = mJ

        # State properties
        # Energy is in cm^-1(NIST database definition)
        self.energy = energy
        self.parity = parity
        
    def degeneracy(self):
        return 2 * self.J +1

    def term_symbol(self):
        L_letters = {
            0: "S",
            1: "P",
            2: "D",
            3: "F",
            4: "G",
        }

        multiplicity = int(2 * self.S + 1)
        L_letter = L_letters.get(self.L, "?")

        return f"{multiplicity}{L_letter}_{self.J}"
    
    def spectroscopic_name(self):
        return f"{self.configuration} {self.term_symbol()}"
    
    
        
class Laser:
    def __init__(self, frequency, intensity, linewidth, polarization = None):
        self.frequency = frequency
        self.intensity = intensity          # intensity unit: W/m^2
        self.linewidth = linewidth
        self.polarization = polarization
        
    def electric_field_amplitude(self):
        return np.sqrt(2 * self.intensity / (c * epsilon_0))
             

class Transition:
    def __init__(self, lower, upper, linewidth=None, lifetime=None, A=None):
        self.lower = lower
        self.upper = upper
        self.linewidth = linewidth
        self.lifetime = lifetime
        self.A = A
        
    def frequency(self):
        return (self.upper.energy - self.lower.energy) * c *100
        
    def rabi_frequency(self, laser):
       E0 = laser.electric_field_amplitude()
       
       return self.dipole * E0 / hbar
    
    def B21(self):
        
        # Need to double check the units
        return self.A * c**3 / (8 * np.pi * h * self.frequency()** 3)
   
    def B12(self):
        g1 = self.lower.degeneracy()
        g2 = self.upper.degeneracy()
        return (g2 / g1) * self.B21()
    
    def detuning(self, laser):
        return self.frequency() - laser.frequency
    
    def transitionDipoleMoment(self):
        d_ij = np.sqrt(
            3 * np.pi * epsilon_0 * hbar * c**3 * self.A
            / self.frequency()**3
            )
        return d_ij
    
        
    
    '''
    # This part is replaced by calculating the Voigt profile directly
    def lorentzian(self, laser):
        return 1 / (1 + (2 * self.detuning(laser) / self.linewidth)**2)
    
    def gaussian(self, laser, sigma):
        return np.exp(-(self.detuning(laser)**2) / (2 * sigma**2))
    '''
    
    def doppler_sigma(self, temperature, mass):
        return (self.frequency() / c * np.sqrt(k * temperature / mass))
    
    
    def line_profile(self, laser, temperature, mass):
         delta = self.detuning(laser)
         sigma = self.doppler_sigma(temperature, mass)
         gamma = self.linewidth / (2 * np.pi)
         
         V = voigt_profile(delta, sigma, gamma)
         V0 = voigt_profile(0.0, sigma, gamma)

         return V / V0
    
class RateEquationSolver:
    def __init__(self, states, transitions, lasers, temperature, mass):
        self.states = states
        self.transitions = transitions
        self.lasers = lasers
        self.temperature = temperature
        self.mass = mass
        
        # s is used as the key in the dictionary, and i is the matrix index
        self.index = {id(s): i for i, s in enumerate(states)}
        
    
    def build_2_step_matrix(self):
        transition12 = self.transitions[0]
    
        laser1 = self.lasers[0]
        laser2 = self.lasers[1]
        
        V12 = transition12.line_profile(laser1, self.temperature, self.mass)
        
        A21 = transition12.A
        
        row12 = laser1.intensity / (c * laser1.linewidth)
        b12v12 = transition12.B12() * V12 * row12
        b21v12 = transition12.B21() * V12 * row12
        
        sigma_ion = 6e-21
        W_ion = sigma_ion * laser2.intensity / (h * laser2.frequency)
        
        M_2_step = np.array([
                    [-b12v12, b21v12 + A21, 0],
                    [-A21 + b12v12, -b21v12 - W_ion, 0 ],
                    [0, W_ion, 0]])
        
        
        return M_2_step
    
    def build_3_step_matrix(self):
        transition12 = self.transitions[0]
        transition23 = self.transitions[1]
        transition13 = self.transitions[2]
    
        laser1 = self.lasers[0]
        laser2 = self.lasers[1]
        laser3 = self.lasers[2]  
    
        '''
        V12 = transition12.line_profile(laser1, self.temperature, self.mass) \
          * (laser1.intensity / c) * 2 / (np.pi * transition12.linewidth)
        V23 = transition23.line_profile(laser2, self.temperature, self.mass) \
          * (laser2.intensity / c) * 2 / (np.pi * transition23.linewidth)
        '''
        
        V12 = transition12.line_profile(laser1, self.temperature, self.mass)
        V23 = transition23.line_profile(laser2, self.temperature, self.mass)
        print(V12)
        print(V23)
        A21 = transition12.A
        A32 = transition23.A
        A31 = transition13.A
        
        row12 = laser1.intensity / (c * laser1.linewidth)
        row23 = laser2.intensity / (c * laser2.linewidth)
        b12v12 = transition12.B12() * V12 * row12
        b21v12 = transition12.B21() * V12 * row12
        b23v23 = transition23.B12() * V23 * row23
        b32v23 = transition23.B21() * V23 * row23
        
        # Need to check later!
        sigma_ion = 6e-21
        W_ion = sigma_ion * laser3.intensity / (h * laser3.frequency)
        print(W_ion)
        
        M = np.array([
                     [-b12v12,                   b21v12 + A21,              A31,    0],
                     [ b12v12,  -(b21v12 + b23v23 + A21),         b32v23 + A32,    0],
                     [      0,                   b23v23,   -(b32v23 + A32 + A31 + W_ion),  0],
                     [      0,                        0,                  W_ion,    0]
        ])
        return M

    '''
    will work on the extensible method latter...
    def build_matrix(self):
        n = len(self.states) + 1
        M = np.zeros((n, n))

        for transition in self.transitions:
            i = self.index[id(transition.lower)]
            j = self.index[id(transition.upper)]

            # Add A coefficient terms
            if transition.A21 is not None:
                A = transition.A21
                M[i, j] += A    
                M[j, j] -= A
            
        # Add B coefficient terms
        for laser in self.lasers:
            V = transition.line_profile(laser, self.temperature, self.mass)
            print(V)
            rate_abs = transition.B12() * V
            rate_emi = transition.B21() * V

            M[j, i] += rate_abs
            M[i, i] -= rate_abs
            M[i, j] += rate_emi   
            M[j, j] -= rate_emi
        print(M)
        return M
    '''
    
    def population_at_t(self, N0, t):
        M = self.build_3_step_matrix()
        return np.ndarray.tolist(expm(M * t) @ N0)
        
    def solve(self, N0, t):
        M = self.build_3_step_matrix()
        return np.array([expm(M * ti) @ N0 for ti in t])
    
    def population_at_t_2step(self, N0, t):
        M = self.build_2_step_matrix()
        return np.ndarray.tolist(expm(M * t) @ N0)
    
    def solve_2step(self, N0, t):
        M = self.build_2_step_matrix()
        print(M)
        return np.array([expm(M * ti) @ N0 for ti in t])

"""
# Gonna use QuTip for this!
class densityMatrixSolver:
    def __init__(self, states, transitions, lasers, mass):
        self.states = states
        self.transitions = transitions
        self.lasers = lasers
        self.mass = mass
        
    def buildHamiltonian(self):
        
    
    
    def buildDissipation(self, densityMatrix):
        
        
        
    def solve(self, densityMatrix, t):
        
"""
    
    
            
            
        
    
    