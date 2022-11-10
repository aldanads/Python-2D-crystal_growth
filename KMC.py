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
    coord_Mo = MoS2_lattice.coord_xy_defects # Coordinates of Mo adatoms
    atomic_specie = MoS2_lattice.atomic_specie # Defect species: 3 (Mo) or 1 (S)
    defect_specie = MoS2_lattice.defect_specie
    
    for i in np.arange(len(coord_Mo[0])):

    
        Mo_adatom = Defects(coord_Mo[0][i],coord_Mo[1][i],MoS2_lattice.Act_E,atomic_specie)
        MoS2_crystal = Cluster(Grid_states)

        TR =Mo_adatom.TR(T,Grid_states,MoS2_crystal.join_cluster_ij) # Calculate the transition rates
        
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

