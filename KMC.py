# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:07:32 2022

@author: ALDANADS
"""
import numpy as np
from defects import*


def KMC(MoS2_lattice, MoS2_crystal,distribution_parameters,prob):
    
    Grid_states = MoS2_lattice.Grid_states # 
    T = MoS2_lattice.T # Crystal temperature 
    
    """
    ----- Coordinates and labels of adatoms and defects at the crystal edge ---
    """
    coord_Mo = MoS2_lattice.coord_defects() # Search the coordinates of every defects - Update at every time step
    cluster_ij = MoS2_crystal.cluster_ij # Edge of the cluster
    xy_adatom_edge = coord_Mo+cluster_ij # Edge of the cluster and adatoms
    l_defect_species = np.zeros(len(xy_adatom_edge))+2
    l_defect_species[len(coord_Mo):] = 4 # Label of defects (adatoms = 2 and edge of the cluster = 4)
    """
    ---------------------------------------------------------------------------
    """
    
    tmax = 0
    pair_atom_defect = distribution_parameters[3]


    if MoS2_lattice.n_defects < 10000:
        distribution = distribution_parameters[0]
        skewness = distribution_parameters[1] 
        fissure_region = distribution_parameters[2] 
        prob_defects = distribution_parameters[4]
  

        MoS2_lattice.defect_distributions(prob_defects,fissure_region,skewness,distribution,pair_atom_defect)


    for i in np.arange(len(xy_adatom_edge)):

    
        Mo_adatom = Defects(xy_adatom_edge[i][0],xy_adatom_edge[i][1],MoS2_lattice.Act_E,pair_atom_defect[0],T,Grid_states,MoS2_crystal.join_cluster_ij,l_defect_species[i])
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
            
        #Calculate the time
        time =  -np.log(np.random.rand())/sum_TR
        if time > tmax:
            tmax=time
            
        # The probability of an event happening at a specific time
        event_prob=1.0-np.exp(-TR[s]*tmax);

        """
        --------------------------------------------------------------------
        """
                
        if np.random.rand() < event_prob:
            # Update the Grid with the chosen events
            Grid_states = Mo_adatom.processes(MoS2_crystal.Grid_states,s) 
            prob[s+1] += 1
            prob[0] += 1

            if (s+1 >= 7): # Update the crystal -> New adatom joined
                Grid_states = MoS2_crystal.crystal_update(Grid_states,(xy_adatom_edge[i][0],xy_adatom_edge[i][1]),s,(Mo_adatom.i,Mo_adatom.j))
                                
        #MoS2_lattice.plot_lattice()

        
        
    MoS2_lattice.Grid_states = Grid_states # Store the new lattice state
    MoS2_lattice.add_time(tmax)
    #print(prob/prob[0])
    return MoS2_lattice,MoS2_crystal,Mo_adatom,prob

