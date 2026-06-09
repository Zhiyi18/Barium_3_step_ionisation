# -*- coding: utf-8 -*-
"""
Created on Fri May 29 11:07:43 2026

@author: USER-01
"""

import requests
import pandas as pd
from io import StringIO

def query_nist_lines(element, low_nm, high_nm):
    """
    Fetches lines from NIST ASD for a given element and wavelength range.
    Returns a pandas DataFrame with wavelength, Aki (A coefficient), 
    Ei/Ek (energy levels), configurations.
    """
    url = "https://physics.nist.gov/cgi-bin/ASD/lines1.pl"
    params = {
        'spectra': element,        # e.g. 'Ba I', 'Ba II'
        'low_w': low_nm,
        'upp_w': high_nm,
        'unit': 1,                 # wavelength in nm
        'format': 2,               # tab-delimited ASCII
        'en_unit': 0,              # energy in cm⁻¹
        'A_out': 0,                # include A coefficients
        'allowed_out': 1,
        'forbid_out': 1,
        'no_spaces': 'on',
        'tab_delimited': 'on',
    }
    
    print("Sending request...")
    response = requests.get(url, params=params, timeout=30)
    print("Got response")

    # Parse tab-delimited response into DataFrame
    df = pd.read_csv(StringIO(response.text), sep='\t', 
                     comment='#', skip_blank_lines=True)
    return df

# Example: get all Ba I lines around 553 nm
df = query_nist_lines('Ba I', 500, 600)
print(df.columns)   # see what columns come back
print(df.head())