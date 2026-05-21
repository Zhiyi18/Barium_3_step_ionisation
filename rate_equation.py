# -*- coding: utf-8 -*-
"""
Created on Thu May 21 11:49:00 2026

@author: USER-01
"""

import numpy as np
import library

ground_state = library.State(
    configuration="6s6p",
    S=0,
    L=1,
    J=1,
    energy=18060
)

print(ground_state.spectroscopic_name())