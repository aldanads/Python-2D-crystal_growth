# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:07:32 2022

@author: ALDANADS
"""
import numpy as np
from defects import Defects
import sys
from balanced_tree import Node, build_tree, update_data, search_value


def KMC(MoS2_lattice, MoS2_crystal,distribution_parameters,events,rng):
    
    Grid_states = MoS2_lattice.Grid_states # 
    T = MoS2_lattice.T # Crystal temperature 

    
    """
    ----- Coordinates and labels of adatoms and defects at the crystal edge ---
    """
    coord_Mo = MoS2_lattice.coord_defects() # Search the coordinates of every defects - Update at every time step
    edge_ij = MoS2_crystal.edge_ij # Atoms at edge of the cluster
    xy_adatom_edge = coord_Mo+edge_ij # Atoms at edge of the cluster and adatoms
    l_defect_species = np.zeros(len(xy_adatom_edge))+2
    l_defect_species[len(coord_Mo):] = 4 # Label of defects (adatoms = 2 and atoms at edge of the cluster = 4)
    """
    ---------------------------------------------------------------------------
    """
    
    time = 0
    previous_step_desorption = events[1][8]
    pair_atom_defect = distribution_parameters[3]
    
    MoS2_lattice.defect_distributions(distribution_parameters)
    
    # Calculate the transition rates of all the events
    TR = []
    for i in np.arange(len(xy_adatom_edge)):

        Mo_adatom = Defects(xy_adatom_edge[i][0],xy_adatom_edge[i][1],MoS2_lattice.Backup_energy,pair_atom_defect[0],T,Grid_states,MoS2_crystal.join_cluster_ij,l_defect_species[i])
        # =============================================================================
        #         TR[0] = Transition rate
        #         TR[1] = Type of event
        #         TR[2] = Particle selected
        # =============================================================================
        TR.extend([(Mo_adatom.TR[j],j, i) for j in np.arange(len(Mo_adatom.TR)) if Mo_adatom.TR[j] != 0.0])

    """
    ------------- Select event with kinetic Monte Carlo technique ------
    """
    # Sort the list of events
    sorted(TR,key = lambda x:x[0])
    # Build a balanced tree structure
    TR_tree = build_tree(TR)
    # Each node is the sum of their children, starting from the leaf
    sumTR = update_data(TR_tree)
    
    # We search in our binary tree the events that happen
    chosen_event = search_value(TR_tree,sumTR*rng.random())

    #Calculate the time step
    time =  -np.log(rng.random())/sumTR

    # events[0][1:] += Mo_adatom.allowed_events[1:]

    # Update the Grid with the chosen events
    Grid_states = Mo_adatom.processes(MoS2_crystal.Grid_states,chosen_event[1],xy_adatom_edge[chosen_event[2]][0],xy_adatom_edge[chosen_event[2]][1],l_defect_species[chosen_event[2]]) 
    # events[1][s+1] += 1
    # events[1][0] += 1

    print(chosen_event[1]+1)
    MoS2_lattice.Grid_states = Grid_states # Store the new lattice state
    MoS2_lattice.coord_defects()
    MoS2_lattice.plot_lattice()
    sys.exit("Error message")

    if (chosen_event[1]+1 >= 7): # Update the crystal -> New adatom joined
        Grid_states = MoS2_crystal.crystal_update(Grid_states,(xy_adatom_edge[i][0],xy_adatom_edge[i][1]),chosen_event[1],(Mo_adatom.i,Mo_adatom.j))
 
        
    #MoS2_lattice.plot_lattice()
        
        
        
    MoS2_lattice.Grid_states = Grid_states # Store the new lattice state
    MoS2_lattice.add_time(time)
    MoS2_lattice.net_flux(events[1][8]-previous_step_desorption)
    MoS2_crystal.crystal_area() # Register the crystal size at this time
    return MoS2_lattice,MoS2_crystal,Mo_adatom,events

