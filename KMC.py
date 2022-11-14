# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:07:32 2022

@author: ALDANADS
"""
import numpy as np
from defects import*

def KMC(MoS2_lattice, MoS2_crystal,pair_atom_defect):
    
    Grid_states = MoS2_lattice.Grid_states # 
    T = MoS2_lattice.T # Crystal temperature 
    coord_Mo = MoS2_lattice.coord_defects() # Search the coordinates of every defects - Update at every time step

    
    for i in np.arange(len(coord_Mo)):

    
        Mo_adatom = Defects(coord_Mo[i][0],coord_Mo[i][1],MoS2_lattice.Act_E,pair_atom_defect[0],T,Grid_states,MoS2_crystal.join_cluster_ij)
        TR = Mo_adatom.TR # Transition rates
        
        """
        ------------- Select event with kinetic Monte Carlo technique ------
        """
        sum_TR = sum(TR)*np.random.rand()
        if sum_TR == 0: continue
        Pointer_event = TR[0]
        s = 0

        while (Pointer_event <= sum_TR):
            s += 1
            Pointer_event += TR[s]
        """
        --------------------------------------------------------------------
        """
        
        # Update the Grid with the chosen events
        Grid_states = Mo_adatom.processes(MoS2_crystal.Grid_states,s) 
        
        
        if (s+1) >= 7: # Update the crystal -> New adatom joined
            Grid_states = MoS2_crystal.crystal_growth(Grid_states,(coord_Mo[i][0],coord_Mo[i][1]))

    MoS2_lattice.Grid_states = Grid_states # Store the new lattice state

    return MoS2_lattice,MoS2_crystal,Mo_adatom

