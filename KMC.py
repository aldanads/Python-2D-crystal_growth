# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:07:32 2022

@author: ALDANADS
"""
import numpy as np


def KMC(Grid_states):
    
    return Grid_states

def mig_available(Grid_states):
    
    return

def TR(Ea,T):
    
    kb = 8.6173324E-5 # Boltzmann constant
    nu0=7E13;  # nu0 (s^-1) bond vibration frequency
    
    return nu0*np.exp(-Ea/(kb*T));