# -*- coding: utf-8 -*-
"""
Spyderエディタ

これは一時的なスクリプトファイルです。
"""

import numpy as np
import matplotlib as plt
from scipy.constants import h, hbar, epsilon_0, c, k

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
        self.energy = energy
        self.parity = parity
        
    def degeneracy(self):
        return 2J +1

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
    def __init__(self, frequency, intensity, polarization):
        self.frequency = frequency
        self.intensity = intensity          # intensity unit: W/m^2
        self.polariztion = polarization
        
    def electric_field_amplitude(self):
        return np.sqrt(2 * self.intensity / (c * epsilon_0))
             

class Transition:
    def __init__(self, lower, upper, linewidth=None, dipole=None, A21=None):
        self.lower = lower
        self.upper = upper
        self.linewidth = linewidth
        self.dipole = dipole
        self.A21 = A21
        
    def frequency(self):
        return (self.upper.energy - self.lower.energy) / h
        
    def rabi_frequency(self, laser):
       E0 = laser.electric_field_amplitude()
       
       return self.dipole * E0 / hbar
    
    def B21(self, laser):
        freq = laser.frequency
        
        # Need to double check the units
        return self.A21 * c**3 / (8 * pi * h * freq**3)
   
    # Is (self, laser) needed here?
    def B12(self, laser):
        g1 = self.lower.degeneracy
        g2 = self.upper.degeneracy
        return (g2 / g1) * self.B21(laser)
    
    def detuning(self, laser):
        return self.frequency - laser.frequency



def lorentzian(detuning, linewidth):
    return 1 / (1 + (2 * detuning / linewidth)**2)

def gaussian(detuning, sigma):
    return np.exp(-(detuning**2) / (2 * sigma**2))

def doppler_sigma(transition, temperature, mass):
    return (transition.frequency / c * np.sqrt(k * temperature / mass))


def line_profile(laser, transition, temperature, mass, frequency):

    detuning = frequency - transition.frequency

    # Check if this approximation is correct!
    gamma = transition.linewidth + laser.linewidth

    sigma = doppler_sigma(transition, temperature, mass)

    L = lorentzian(detuning, gamma)
    G = gaussian(detuning, sigma)

    return L * G
        
    
def rate_equation_solver(line1, line2,
                         laser1, laser2, laser3):
    
    