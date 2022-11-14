# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:07:32 2022

@author: ALDANADS
"""
import numpy as np
from defects import*

def KMC(MoS2_lattice, MoS2_crystal):
    
    Grid_states = MoS2_lattice.Grid_states # 
    T = MoS2_lattice.T # Crystal temperature 
    coord_Mo = MoS2_lattice.coord_defects() # Search the coordinates of every defects - Update at every time step

    
    for i in np.arange(len(coord_Mo)):

    
        Mo_adatom = Defects(coord_Mo[i][0],coord_Mo[i][1],MoS2_lattice.Act_E,MoS2_lattice.atomic_specie)
        #"""
        MoS2_crystal = Cluster(Grid_states) #--> This should be an update
        #"""
        
        TR = Mo_adatom.TR(T,Grid_states,MoS2_crystal.join_cluster_ij) # Calculate the transition rates
        
        # Select event
        sum_TR = sum(TR)*np.random.rand()
        if sum_TR == 0: continue
        Pointer_event = TR[0]
        s = 0

        while (Pointer_event <= sum_TR):
            s += 1
            Pointer_event += TR[s]
        
        
        Grid_states = Mo_adatom.processes(MoS2_crystal.Grid_states,s) # Update the Grid with new events

    MoS2_lattice.Grid_states = Grid_states # Store the new lattice state

    return MoS2_lattice,MoS2_crystal

