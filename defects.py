# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 14:03:55 2022

@author: ALDANADS
"""
import numpy as np

class Defects():
    
    
    def __init__(self,i,j,Act_E,atomic_specie):
        
        self.i = i
        self.j = j
        self.Act_E = Act_E
        self.atomic_specie = atomic_specie

    
    """
    ________________________ EVENTS AVAILABLE __________________________
    """
    
    # Check possible events
    def events_available(self,Grid_states):
        

        
        # Kink nucleation
        
        # Kink propagation
        
        """ Adatom migration: Mo migrating to Mo positions
            # Mo left (2) or Mo right (3) allowed_events[0]
            # Down - Zigzag - allowed_events[1]    
            # Up - Zigzag - allowed_events[2]   
            # Left up - Armchair - allowed_events[3]
            # Left down - Armchair - allowed_events[4]
            # Right up - Armchair - allowed_events[5]
            # Right down - Armchair - allowed_events[6]

        """
        allowed_events = np.zeros(7)
        self.allowed_events = allowed_events 
        # allowed_events =
        # 0: forbidden
        # 1: allowed
        # 2: left Mo
        # 3: right Mo
        
        self.mig_available(Grid_states) # Migration
        
        # Function: Check possible migration routes
    def mig_available(self,Grid_states):
        
        length_x = len(Grid_states)-1
        length_y = len(Grid_states[0])-1
        
        i = self.i
        j = self.j
        
        atomic_specie = self.atomic_specie
        
        if ((i>0) and (Grid_states[i-1,j] == 0)) or ((i<length_x) and (Grid_states[i+1,j] == 0)):
        # 0
        # 2: Mo (Grid_states) Left Mo --------------------------------------------
        # 0
        
            self.allowed_events[0] = 2 # Left Mo
            """
            Up and down
            """
            if (i>1) and Grid_states[i-2,j] == atomic_specie: # Down
               self.allowed_events[1] = 1 
               
            if (i<length_x-1) and (Grid_states[i+2,j] == atomic_specie): # Up
                self.allowed_events[2] = 1
                
            """
            Left up and down
            """
            if (j>1) and (i<length_x) and (Grid_states[i+1,j-2] == atomic_specie): # Left up
                self.allowed_events[3] = 1

            if (j>1) and (i>0) and (Grid_states[i-1,j-2] == atomic_specie): # Left down
                self.allowed_events[4] = 1
                
            """
            Right up and down
            """
            if (i<length_x) and (j<length_y) and (Grid_states[i+1,j+1] == atomic_specie): # Right up
                self.allowed_events[5] = 1
                
            if (i>0) and (j<length_y) and (Grid_states[i-1,j+1] == atomic_specie): # Right down
                self.allowed_events[6] = 1

                
        
        if ((i>0) and (Grid_states[i-1,j] == 1)) or ((i<length_x) and (Grid_states[i+1,j] == 1)):

        # S
        # 3: Mo (Grid_states) Right Mo --------------------------------------------
        # S
            self.allowed_events[0] = 3 # Right Mo
            
            """
            Up and down
            """
            if (i>1) and Grid_states[i-2,j] == atomic_specie: # Down
               self.allowed_events[1] = 1 
               
            if (i<length_x-1) and (Grid_states[i+2,j] == atomic_specie): # Up
                self.allowed_events[2] = 1

            """
            Left up and down
            """
            if (j>0) and (i<length_x) and (Grid_states[i+1,j-1] == atomic_specie): # Left up
                self.allowed_events[3] = 1
            
            if (j>0) and (i>0) and (Grid_states[i-1,j-1] == atomic_specie): # Left down
                self.allowed_events[4] = 1
                
            """
            Right up and down
            """
            if (i<length_x) and (j<length_y-1) and (Grid_states[i+1,j+2] == atomic_specie): # Right up
                self.allowed_events[5] = 1
                    
            if (i>0) and (j<length_y-1) and (Grid_states[i-1,j+2] == atomic_specie): # Right down
                self.allowed_events[6] = 1
        
        
        """
        ________________________ TRANSITION RATES __________________________
        """
        
        
    def TR(self,T,Grid_states):
        
        self.events_available(Grid_states)
        kb = 8.6173324E-5 # Boltzmann constant
        nu0=7E13;  # nu0 (s^-1) bond vibration frequency
        allowed_events = self.allowed_events
        TR = np.zeros(len(allowed_events)-1)
            
        TR = nu0*np.exp(-np.array(self.Act_E)/(kb*T))
        TR[allowed_events[1:] == 0] = 0
        self.TR = TR
        
        """
        ________________________ PROCESSES __________________________
        """
        

    def processes(self,Grid_states,s):
        
        s = s+1 # s is selected from TR, which is smaller than allowed_events
        
        i = self.i
        j = self.j
        
        atomic_specie = self.atomic_specie
        
        # Down - Zigzag - allowed_events[1]    
        # Up - Zigzag - allowed_events[2]   
        # Left up - Armchair - allowed_events[3]
        # Left down - Armchair - allowed_events[4]
        # Right up - Armchair - allowed_events[5]
        # Right down - Armchair - allowed_events[6]
        
        if self.allowed_events[0] == 2: # Mo left (2) 
        
            """
            ---------------------- Migration ----------------------------------
            """
            
        
            if (s == 1): # Down - Zigzag
                Grid_states[i,j] = atomic_specie
                Grid_states[i-2,j] = 2
                self.i = i-2
                
            elif (s == 2): # Up - Zigzag
                Grid_states[i,j] = atomic_specie
                Grid_states[i+2,j] = 2
                self.i = i+2
                
            elif (s == 3): # Left up - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i+1,j-2] = 2
                self.i = i+1
                self.j = j-2
            
            elif (s == 4): # Left down - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i-1,j-2] = 2
                self.i = i-1
                self.j = j-2
                
            elif (s == 5): # Right up - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i+1,j+1] = 2
                self.i = i+1
                self.j = j+1
                
            elif (s == 6): # Right down - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i-1,j+1] = 2
                self.i = i-1
                self.j = j+1
        
        elif self.allowed_events[0] == 3: # Mo right (3)
        
            """
            ---------------------- Migration ----------------------------------
            """
            
            if (s == 1): # Down - Zigzag
                Grid_states[i,j] = atomic_specie
                Grid_states[i-2,j] = 2
                self.i = i-2
                
            elif (s == 2): # Up - Zigzag
                Grid_states[i,j] = atomic_specie
                Grid_states[i+2,j] = 2
                self.i = i+2
                
            elif (s == 3): # Left up - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i+1,j-1] = 2
                self.i = i+1
                self.j = j-1
            
            elif (s == 4): # Left down - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i-1,j-1] = 2
                self.i = i-1
                self.j = j-1
            
            elif (s == 5): # Right up - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i+1,j+2] = 2
                self.i = i+1
                self.j = j+2
            
            elif (s == 6): # Right down - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i-1,j+2] = 2
                self.i = i-1
                self.j = j+2
            
        return Grid_states
    




            
class Cluster():
        
    def __init__(self,Grid_states,atomic_specie):
            
        self.cluster_ij = np.where(Grid_states == 4)
        self.cluster_size = sum(sum(Grid_states == 4)) + sum(sum(Grid_states == 5))
        self.atomic_specie = atomic_specie
        self.clustering_region(Grid_states)
            
        
    # Search for the region where the adatom can join the growing crystal
    def clustering_region(self,Grid_states):
            

            
        length_x = len(Grid_states)-1
        length_y = len(Grid_states[0])-1
                    
        atomic_specie = self.atomic_specie
                    
        # Create a empty tuple (with two lists) for the region where the 
        # joining to the cluster is possible
        join_cluster_ij = []   
        
        # For loop over the all the particles in the crystal edge
        for k in np.arange(len(self.cluster_ij[0])):    
            i = self.cluster_ij[0][k]
            j = self.cluster_ij[1][k]
            join_sites = 0 # Number of free sites around a specific cluster point
                
            """
            //////////////////////////////// Left Mo ////////////////////////
            """
            if ((i>0) and (Grid_states[i-1,j] == 0)) or ((i<length_x) and (Grid_states[i+1,j] == 0)):
    
                """
                Up and down
                """
                if (i>1) and (Grid_states[i-2,j] == atomic_specie): # Down
                    join_cluster_ij.append((i-2,j))    
                    join_sites += 1
                        
                if (i<length_x-1) and (Grid_states[i+2,j] == atomic_specie): # Up
                    join_cluster_ij.append((i+2,j))    
                    join_sites += 1
                        
                """
                Left up and down
                """
                if (j>1) and (i<length_x) and (Grid_states[i+1,j-2] == atomic_specie): # Left up
                    join_cluster_ij.append((i+1,j-2))    
                    join_sites += 1
                            
                if (j>1) and (i>0) and (Grid_states[i-1,j-2] == atomic_specie): # Left down
                    join_cluster_ij.append((i-1,j-2))    
                    join_sites += 1
                        
                """
                Right up and down
                """
                if (i<length_x) and (j<length_y) and (Grid_states[i+1,j+1] == atomic_specie): # Right up
                    join_cluster_ij.append((i+1,j+1))    
                    join_sites += 1
                                
                if (i>0) and (j<length_y) and (Grid_states[i-1,j+1] == atomic_specie): # Right down
                    join_cluster_ij.append((i-1,j+1))    
                    join_sites += 1
                        
                """
                //////////////////////////////// Right Mo ////////////////////////
                """
                
            if ((i>0) and (Grid_states[i-1,j] == 1)) or ((i<length_x) and (Grid_states[i+1,j] == 1)):
    
                """
                Up and down
                """
                if (i>1) and (Grid_states[i-2,j] == atomic_specie): # Down
                    join_cluster_ij.append((i-2,j))    
                    join_sites += 1
                        
                if (i<length_x-1) and (Grid_states[i+2,j] == atomic_specie): # Up
                    join_cluster_ij.append((i+2,j))    
                    join_sites += 1
                    
                """
                Left up and down
                """
                if (j>0) and (i<length_x) and (Grid_states[i+1,j-1] == atomic_specie): # Left up
                    join_cluster_ij.append((i+1,j-1))    
                    join_sites += 1
                        
                if (j>0) and (i>0) and (Grid_states[i-1,j-1] == atomic_specie): # Left down
                    join_cluster_ij.append((i-1,j-1))    
                    join_sites += 1
                        
                """
                Right up and down
                """
                if (i<length_x) and (j<length_y-1) and (Grid_states[i+1,j+2] == atomic_specie): # Right up
                    join_cluster_ij.append((i+1,j+2))    
                    join_sites += 1
                        
                if (i>0) and (j<length_y-1) and (Grid_states[i-1,j+2] == atomic_specie): # Right down
                    join_cluster_ij.append((i-1,j+2))    
                    join_sites += 1
                    
            if join_sites == 0: Grid_states[i,j] = 5 # Inner point of the crystal
            
        # Grid points adjacent to the crystal, so they may join the crystal
        self.join_cluster_ij = set(join_cluster_ij) # Select unique elements from the list
        self.Grid_states = Grid_states
                        

                  
 
