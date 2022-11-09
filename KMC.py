# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:07:32 2022

@author: ALDANADS
"""
import numpy as np
from defects import*

def KMC(MoS2_lattice):
    
    Grid_states = MoS2_lattice.Grid_states # 
    T = MoS2_lattice.T # Crystal temperature 
    MoS2_lattice.coord_defects() # Search the coordinates of every defects
    coord_Mo = MoS2_lattice.coord_xy_defects
    for i in np.arange(len(coord_Mo[0])):

    
        Mo_adatom = Defects(coord_Mo[0][i],coord_Mo[1][i],MoS2_lattice.Act_E)
        #Mo_adatom = Defects(0,0,MoS2_lattice.Act_E)
        
        Mo_adatom.TR(T,Grid_states)
        TR = Mo_adatom.TR
        
        # Select event
        sum_TR = sum(TR)*np.random.rand()
        if sum_TR == 0: continue
        Pointer_event = TR[0]
        s = 0
        
        while (Pointer_event <= sum_TR):
            s += 1
            Pointer_event += TR[s]
            
        #print(s+1)
        
        Grid_states = Mo_adatom.processes(Grid_states,s)
        
    MoS2_lattice.Grid_states = Grid_states
    MoS2_lattice.plot_lattice()

    return MoS2_lattice, Mo_adatom

