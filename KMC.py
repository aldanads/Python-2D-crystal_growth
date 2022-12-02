# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:07:32 2022

@author: ALDANADS
"""
import numpy as np
from defects import Defects


def KMC(MoS2_lattice, MoS2_crystal,distribution_parameters,events,rng):
    
    Grid_states = MoS2_lattice.Grid_states # 
    T = MoS2_lattice.T # Crystal temperature 

    
    """
    ----- Coordinates and labels of adatoms and defects at the crystal edge ---
    """
    coord_Mo = MoS2_lattice.coord_defects() # Search the coordinates of every defects - Update at every time step
    edge_ij = MoS2_crystal.edge_ij # Edge of the cluster
    xy_adatom_edge = coord_Mo+edge_ij # Edge of the cluster and adatoms
    l_defect_species = np.zeros(len(xy_adatom_edge))+2
    l_defect_species[len(coord_Mo):] = 4 # Label of defects (adatoms = 2 and edge of the cluster = 4)
    """
    ---------------------------------------------------------------------------
    """
    
    time = 0
    pair_atom_defect = distribution_parameters[3]
    
    MoS2_lattice.defect_distributions(distribution_parameters)


    for i in np.arange(len(xy_adatom_edge)):

    
        Mo_adatom = Defects(xy_adatom_edge[i][0],xy_adatom_edge[i][1],MoS2_lattice.Backup_energy,pair_atom_defect[0],T,Grid_states,MoS2_crystal.join_cluster_ij,l_defect_species[i])
        TR = Mo_adatom.TR # Transition rates
        
        """
        ------------- Select event with kinetic Monte Carlo technique ------
        """
        sum_TR = sum(TR)*rng.random()
        

        
        if sum_TR == 0: continue
        pointer_event = TR[0]
        s = 0

        while (pointer_event < sum_TR):
            s += 1
            pointer_event += TR[s]


        #Calculate the time
        time =  -np.log(rng.random())/sum(TR)
        # if time > tmax:
        #     tmax=time
            
        # The probability of an event happening at a specific time
        event_prob=1.0-np.exp(-TR[s]*time);
        

        """
        --------------------------------------------------------------------
        """
        events[0][1:] += Mo_adatom.allowed_events[1:]
        
        if rng.random() < event_prob:
            

            # Update the Grid with the chosen events
            Grid_states = Mo_adatom.processes(MoS2_crystal.Grid_states,s) 
            events[1][s+1] += 1
            events[1][0] += 1


            if (s+1 >= 7): # Update the crystal -> New adatom joined
                Grid_states = MoS2_crystal.crystal_update(Grid_states,(xy_adatom_edge[i][0],xy_adatom_edge[i][1]),s,(Mo_adatom.i,Mo_adatom.j))
     
                
        #MoS2_lattice.plot_lattice()
        
        
        
    MoS2_lattice.Grid_states = Grid_states # Store the new lattice state
    MoS2_lattice.add_time(time)
    MoS2_crystal.crystal_area() # Register the crystal size at this time
    return MoS2_lattice,MoS2_crystal,Mo_adatom,events

