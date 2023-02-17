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
    
    
        
    time = 0
    previous_step_desorption = events[8]
    pair_atom_defect = distribution_parameters[3]
    
    MoS2_lattice.defect_distributions(distribution_parameters)
    
    Grid_states = MoS2_lattice.Grid_states # 
    T = MoS2_lattice.T # Crystal temperature 

    
    """
    ----- Coordinates and labels of adatoms and defects at the crystal edge ---
    """
    coord_Mo = MoS2_lattice.coord_defects() # Search the coordinates of every defects - Update at every time step
    edge_ij = MoS2_crystal.edge_ij # Atoms at edge of the cluster
    xy_adatom_edge = coord_Mo+edge_ij # Atoms at edge of the cluster and adatoms
    """
    ---------------------------------------------------------------------------
    """
    Dict_defects = {(xy_adatom_edge[i][0],xy_adatom_edge[i][1]): None for i in np.arange(len(xy_adatom_edge))}
    
    """
    ------------- Calculate the transition rates of all the events ------
    """
    for k in np.arange(len(Dict_defects)):
        TR = []
        
        # =============================================================================
        # =============================================================================
        #         TR[0] = Transition rate
        #         TR[1] = Type of event
        #         TR[2] = Particle selected
        # =============================================================================
        # =============================================================================
        for key in Dict_defects:
            if Dict_defects[key] == None:
                
                Dict_defects[key] = Defects(key[0],key[1],MoS2_lattice.Backup_energy,pair_atom_defect[0],T,Grid_states,MoS2_crystal.join_cluster_ij)
            
            TR.extend([(Dict_defects[key].TR[j],j,key) for j in np.arange(len(Dict_defects[key].TR)) if Dict_defects[key].TR[j] != 0.0])


        """
        ------------- Select event with kinetic Monte Carlo technique ------
        """
        # Sort the list of events
        sorted(TR,key = lambda x:x[0])
        # Build a balanced tree structure
        TR_tree = build_tree(TR)
        # Each node is the sum of their children, starting from the leaf
        sumTR = update_data(TR_tree)
        
        if sumTR == None: break
        # When we only have one node in the tree, it returns a tuple
        if type(sumTR) is tuple: sumTR = sumTR[0]
        # We search in our binary tree the events that happen
        chosen_event = search_value(TR_tree,sumTR*rng.random())
        #Calculate the time step
        time += -np.log(rng.random())/sumTR
    

        # Update the Grid with the chosen events
        Mo_adatom = Dict_defects[chosen_event[2]]

        # We made a list of particles affected by the event to recalculate those TR
        List_neighbors = Mo_adatom.List_neighbors
        Grid_states = Mo_adatom.processes(MoS2_crystal.Grid_states,chosen_event[1],chosen_event[2][0],chosen_event[2][1]) 
        events[chosen_event[1]+1] += 1
        events[0] += 1

        """
        #Updates --> Can I put all the operations inside the objects?
        """
        """
        ------------- Updates depending on the event ------
        """
        # Crystal update or desorption
        if (chosen_event[1]+1 >= 7): # Update the crystal -> New adatom joined
            Grid_states = MoS2_crystal.crystal_update(Grid_states,(chosen_event[2][0],chosen_event[2][1]),chosen_event[1],(Mo_adatom.i,Mo_adatom.j))
            
            if (chosen_event[1]+1 == 7):
                List_neighbors += (Mo_adatom.neighbors(Mo_adatom.i,Mo_adatom.j,Grid_states) + [(Mo_adatom.i,Mo_adatom.j)])

            #Desorption               
            elif (chosen_event[1]+1 == 8):
                del Dict_defects[chosen_event[2]]
                MoS2_lattice.coord_xy_defects[0].remove(chosen_event[2][0])
                MoS2_lattice.coord_xy_defects[1].remove(chosen_event[2][1])
            # Update of migration of atoms at the edge
            elif (chosen_event[1]+1 >= 9):
                if Grid_states[(Mo_adatom.i, Mo_adatom.j)] == 4:
                    Dict_defects[(Mo_adatom.i, Mo_adatom.j)] = Dict_defects[chosen_event[2]]
                del Dict_defects[chosen_event[2]]
                List_neighbors += (Mo_adatom.neighbors(Mo_adatom.i,Mo_adatom.j,Grid_states) + [(Mo_adatom.i,Mo_adatom.j)])
        else:
            # Adatom migration update
            MoS2_lattice.coord_xy_defects[0].remove(chosen_event[2][0])
            MoS2_lattice.coord_xy_defects[1].remove(chosen_event[2][1])
            MoS2_lattice.coord_xy_defects[0].append(Mo_adatom.i)
            MoS2_lattice.coord_xy_defects[1].append(Mo_adatom.j)
            Dict_defects[(Mo_adatom.i, Mo_adatom.j)] = Dict_defects[chosen_event[2]]
            del Dict_defects[chosen_event[2]]

            List_neighbors += (Mo_adatom.neighbors(Mo_adatom.i,Mo_adatom.j,Grid_states) + [(Mo_adatom.i,Mo_adatom.j)])
            
            
        """
        ------------- We need to recalculate all the TR of particles affected by the event ------
        """
        List_neighbors = list(set(List_neighbors))
        
        for key in List_neighbors:
            Dict_defects[key] = None
            


        MoS2_lattice.Grid_states = Grid_states # Store the new lattice state
        
    
    #sys.exit()

    MoS2_lattice.net_flux(events[8]-previous_step_desorption)
    MoS2_crystal.crystal_area() # Register the crystal size at this time
    MoS2_lattice.add_time(time)
    return MoS2_lattice,MoS2_crystal,events

