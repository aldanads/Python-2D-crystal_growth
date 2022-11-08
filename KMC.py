# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:07:32 2022

@author: ALDANADS
"""
import numpy as np
from defects import*

def KMC(MoS2_lattice):
    
    Grid_states = MoS2_lattice.Grid_states
    MoS2_lattice.coord_defects()
    coord_Mo = MoS2_lattice.coord_xy_defects
    T = MoS2_lattice.T

    Mo_adatom = Defects(coord_Mo[0][0],coord_Mo[1][0],MoS2_lattice.Act_E)
    #Mo_adatom = Defects(0,0,MoS2_lattice.Act_E)
    
    Mo_adatom.TR(T,Grid_states)
    TR = Mo_adatom.TR
    
    # Select event
    sum_TR = sum(TR)*np.random.rand()
    #if sum_TR == 0: continue
    Pointer_event = TR[0]
    s = 0
    
    while (Pointer_event <= sum_TR):
        s += 1
        Pointer_event += TR[s]
    
    return MoS2_lattice, Mo_adatom

