# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:07:32 2022

@author: ALDANADS
"""
import numpy as np


def KMC(Grid_states):
    
    return Grid_states

def events_available(Grid_states):
    
    # Kink nucleation
    
    # Kink propagation
    
    """ Migration:
        # Up - Zigzag
        # Down - Zigzag
        # Right up - Armchair
        # Right down - Armchair
        # Left up - Armchair
        # Left down - Armchair
    """
    
    
    
    return

def TR(Ea,T):
    
    kb = 8.6173324E-5 # Boltzmann constant
    nu0=7E13;  # nu0 (s^-1) bond vibration frequency
    
    return nu0*np.exp(-Ea/(kb*T));