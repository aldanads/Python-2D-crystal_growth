# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:07:32 2022

@author: ALDANADS
"""
import numpy as np
from defects import*


def KMC(MoS2_lattice):
    
    Grid_states = MoS2_lattice.Grid_states
    coord_xy_Vs = np.where(Grid_states == 2)

    #Mo_adatom = Defects(coord_xy_Vs[0][0],coord_xy_Vs[1][0])
    Mo_adatom = Defects(0,30)
    Mo_adatom.events_available(Grid_states)
    return MoS2_lattice, Mo_adatom

def TR(Ea,T):
    
    kb = 8.6173324E-5 # Boltzmann constant
    nu0=7E13;  # nu0 (s^-1) bond vibration frequency
    
    return nu0*np.exp(-Ea/(kb*T));