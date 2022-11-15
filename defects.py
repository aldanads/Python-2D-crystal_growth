# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 14:03:55 2022

@author: ALDANADS
"""
import numpy as np

class Defects():
    
    
    def __init__(self,i,j,Act_E,atomic_specie,T,Grid_states,join_cluster_ij,defect_specie):
        
        self.i = i
        self.j = j
        self.Act_E = Act_E
        self.atomic_specie = atomic_specie
        self.defect_specie = defect_specie
        
        # Calculate the transition rates
        self.TR(T,Grid_states,join_cluster_ij)

    
    """
    ________________________ EVENTS AVAILABLE __________________________
    """
    
    # Check possible events
    def events_available(self,Grid_states,join_cluster_ij):
        
        """ Adatom migration: Mo migrating to Mo positions
            # Mo left (2) or Mo right (3) allowed_events[0]
            # Down - Zigzag - allowed_events[1]    
            # Up - Zigzag - allowed_events[2]   
            # Left up - Armchair - allowed_events[3]
            # Left down - Armchair - allowed_events[4]
            # Right up - Armchair - allowed_events[5]
            # Right down - Armchair - allowed_events[6]
            # Nucleation - allowed_events[7]
            # Propagation - allowed_events[8]
            # Desorption - allowed_events[9]
        """

        allowed_events = np.zeros(len(self.Act_E)+1)

        # Adatoms --> Migrate, nucleate or desorpt
        if self.defect_specie == 2: 
    
            self.allowed_events = allowed_events 
            # allowed_events[0] =
            # 2: left Mo
            # 3: right Mo
            # allowed_events[1:end] =
            # 0: forbidden
            # 1: allowed
            
            self.nucleation_available(join_cluster_ij,Grid_states) # Nucleation/propagation
            
            self.mig_available(Grid_states) # Migration
            
            self.allowed_events[9] = 1
            
        # Atoms at the edge of the crystal --> They can migrate through the edge
        if self.defect_specie == 4:
            
            self.allowed_events = allowed_events 
            self.migration_edge(join_cluster_ij,Grid_states)
            
            
        
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
                        
    def nucleation_available(self,join_cluster_ij,Grid_states):
        
        i = self.i
        j = self.j
        length_x = len(Grid_states)-1
        length_y = len(Grid_states[0])-1

        join_cluster_ij = list(set(join_cluster_ij))
 
        # Defect in (i,j) is in contact with the cluster
        if (i,j) in join_cluster_ij: 
            
            # Kink nucleation --> Armchair
            # Search in a square 3x4 around (i,j) --> Remember that init:final:steps exclude final
            if (i>0 and i < length_x and j > 1 and j < length_y-2) and (Grid_states[i-1:i+2:1,j-2:j+3:1] == 4).any():
                self.allowed_events[7] = 1
            # Kink propagation --> Zigzag
            if ((i>1) and Grid_states[i-2,j]) == 4 or (i<length_x-1) and (Grid_states[i+2,j] == 4):
                self.allowed_events[8] = 1 
                
    def migration_edge(self,join_cluster_ij,Grid_states):

        i = self.i
        j = self.j
        
        # If the atom is surrounded by more than 2 other atoms, it is immobile
        if sum(sum(Grid_states[i-2:i+3:1,j-2:j+3:1] >= 4)) > 4: return

        
        length_x = len(Grid_states)-1
        length_y = len(Grid_states[0])-1
        atomic_specie = self.atomic_specie

    
        # We select the position supported by at least three element of the crystal
        # That is, the defect in the edge can't jump outside the crystal during the migration process
        # Jump outside the crystal is the detach process
        join_cluster_ij = [x for x in join_cluster_ij if join_cluster_ij.count(x) > 2]
        
        
        if ((i>0) and (Grid_states[i-1,j] == 0)) or ((i<length_x) and (Grid_states[i+1,j] == 0)):
        # 0
        # 2: Mo (Grid_states) Left Mo --------------------------------------------
        # 0
                
            self.allowed_events[0] = 2 # Left Mo
            """
            Up and down
            """
            if (i>1) and ((i-2,j) in join_cluster_ij and Grid_states[i-2,j] == atomic_specie): # Down
               self.allowed_events[10] = 1 
               
            if (i<length_x-1) and ((i+2,j) in join_cluster_ij and Grid_states[i+2,j] == atomic_specie): # Up
                self.allowed_events[11] = 1
                
            """
            Left up and down
            """
            if (j>1) and (i<length_x) and ((i+1,j-2) in join_cluster_ij and Grid_states[i+1,j-2] == atomic_specie): # Left up
                self.allowed_events[12] = 1

            if (j>1) and (i>0) and ((i-1,j-2) in join_cluster_ij and Grid_states[i-1,j-2] == atomic_specie): # Left down
                self.allowed_events[13] = 1
                
            """
            Right up and down
            """
            if (i<length_x) and (j<length_y) and ((i+1,j+1) in join_cluster_ij and Grid_states[i+1,j+1] == atomic_specie): # Right up
                self.allowed_events[14] = 1
                
            if (i>0) and (j<length_y) and ((i-1,j+1) in join_cluster_ij and Grid_states[i-1,j+1] == atomic_specie): # Right down
                self.allowed_events[15] = 1
                
                

        if ((i>0) and (Grid_states[i-1,j] == 1)) or ((i<length_x) and (Grid_states[i+1,j] == 1)):
        # S
        # 3: Mo (Grid_states) Right Mo --------------------------------------------
        # S
        
            self.allowed_events[0] = 3 # Right Mo
            
            """
            Up and down
            """
            if (i>1) and ((i-2,j) in join_cluster_ij and Grid_states[i-2,j] == atomic_specie): # Down
               self.allowed_events[10] = 1 
               
            if (i<length_x-1) and ((i+2,j) in join_cluster_ij and Grid_states[i+2,j] == atomic_specie): # Up
                self.allowed_events[11] = 1

            """
            Left up and down
            """
            if (j>0) and (i<length_x) and ((i+1,j-1) in join_cluster_ij and Grid_states[i+1,j-1] == atomic_specie): # Left up
                self.allowed_events[12] = 1
            
            if (j>0) and (i>0) and ((i-1,j-1) in join_cluster_ij and Grid_states[i-1,j-1] == atomic_specie): # Left down
                self.allowed_events[13] = 1
                
            """
            Right up and down
            """
            if (i<length_x) and (j<length_y-1) and ((i+1,j+2) in join_cluster_ij and Grid_states[i+1,j+2] == atomic_specie): # Right up
                self.allowed_events[14] = 1
                    
            if (i>0) and (j<length_y-1) and ((i-1,j+2) in join_cluster_ij and Grid_states[i-1,j+2] == atomic_specie): # Right down
                self.allowed_events[15] = 1
        
        
        """
        ________________________ TRANSITION RATES __________________________
        """
        
        
    def TR(self,T,Grid_states,join_cluster_ij):
        
        self.events_available(Grid_states,join_cluster_ij)
        kb = 8.6173324E-5 # Boltzmann constant
        nu0=7E13;  # nu0 (s^-1) bond vibration frequency
        #nu0 = 1E12 #  Shuai, Chen, Gao Junfeng, Bharathi M. Srinivasan, and Zhang Yong-Wei. "A kinetic Monte Carlo study for mono-and bi-layer growth of MoS2 during chemical vapor deposition." Acta Physico-Chimica Sinica 35, no. 10 (2019): 1119-1127.
        allowed_events = self.allowed_events
        TR = np.zeros(len(allowed_events)-1)
            
        TR = nu0*np.exp(-np.array(self.Act_E)/(kb*T))
        TR[allowed_events[1:] == 0] = 0
        self.TR = TR
        
        return TR
        
        """
        ________________________ PROCESSES __________________________
        """
        

    def processes(self,Grid_states,s):
        
        s = s+1 # s is selected from TR, which is smaller than allowed_events
        
        i = self.i
        j = self.j
        
        atomic_specie = self.atomic_specie
        defect_specie = self.defect_specie
        
        # Down - Zigzag - allowed_events[1]    
        # Up - Zigzag - allowed_events[2]   
        # Left up - Armchair - allowed_events[3]
        # Left down - Armchair - allowed_events[4]
        # Right up - Armchair - allowed_events[5]
        # Right down - Armchair - allowed_events[6]
        # Kink nucleation - allowed_events[7]
        # Kink propagation - allowed_events[8]
        # Desorption - allowed_events[9]
        if (s == 7) or (s == 8): # The defect join the cluster -> Exit the function
            Grid_states[i,j] = 4 
            return Grid_states
        
        if (s == 9):
            Grid_states[i,j] = 3
            return Grid_states
        
        elif self.allowed_events[0] == 2: # Mo left (2) 
        
            """
            ---------------------- Migration ----------------------------------
            """
            
        
            if (s == 1 or s == 10): # Down - Zigzag
                Grid_states[i,j] = atomic_specie
                Grid_states[i-2,j] = defect_specie
                self.i = i-2
                
            elif (s == 2 or s == 11): # Up - Zigzag
                Grid_states[i,j] = atomic_specie
                Grid_states[i+2,j] = defect_specie
                self.i = i+2
                
            elif (s == 3 or s == 12): # Left up - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i+1,j-2] = defect_specie
                self.i = i+1
                self.j = j-2
            
            elif (s == 4 or s == 13): # Left down - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i-1,j-2] = defect_specie
                self.i = i-1
                self.j = j-2
                
            elif (s == 5 or s == 14): # Right up - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i+1,j+1] = defect_specie
                self.i = i+1
                self.j = j+1
                
            elif (s == 6 or s == 15): # Right down - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i-1,j+1] = defect_specie
                self.i = i-1
                self.j = j+1
        
        elif self.allowed_events[0] == 3: # Mo right (3)
        
            """
            ---------------------- Migration ----------------------------------
            """
            
            if (s == 1 or s == 10): # Down - Zigzag
                Grid_states[i,j] = atomic_specie
                Grid_states[i-2,j] = defect_specie
                self.i = i-2
                
            elif (s == 2 or s == 11): # Up - Zigzag
                Grid_states[i,j] = atomic_specie
                Grid_states[i+2,j] = defect_specie
                self.i = i+2
                
            elif (s == 3 or s == 12): # Left up - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i+1,j-1] = defect_specie
                self.i = i+1
                self.j = j-1
            
            elif (s == 4 or s == 13): # Left down - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i-1,j-1] = defect_specie
                self.i = i-1
                self.j = j-1
            
            elif (s == 5 or s == 14): # Right up - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i+1,j+2] = defect_specie
                self.i = i+1
                self.j = j+2
            
            elif (s == 6 or s == 15): # Right down - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i-1,j+2] = defect_specie
                self.i = i-1
                self.j = j+2
            
        return Grid_states
    




            
class Cluster():
        
    def __init__(self,Grid_states):
            
        cluster_ij = np.where(Grid_states == 4)
        cluster_ij = [(cluster_ij[0][i],cluster_ij[1][i]) for i in np.arange(len(cluster_ij[0]))]

        self.cluster_ij = cluster_ij # List of tuples
        self.cluster_size = sum(sum(Grid_states == 4)) + sum(sum(Grid_states == 5))
        
        # Create a empty list --> Add tuples with coordinates where defects can 
        # joining to the cluster
        self.join_cluster_ij = []   
        
        self.clustering_region(Grid_states,cluster_ij)
            
        
    # Search for the region where the adatom can join the growing crystal
    def clustering_region(self,Grid_states,cluster_ij):
            

            
        length_x = len(Grid_states)-1
        length_y = len(Grid_states[0])-1
                    
                    
        join_cluster_ij = self.join_cluster_ij # Need to be a list to append elements
        
        # For loop over the all the particles in the crystal edge
        for k in np.arange(len(cluster_ij)): 

            i = cluster_ij[k][0]
            j = cluster_ij[k][1]
            join_sites = 0 # Number of free sites around a specific cluster point
                
            """
            //////////////////////////////// Left Mo ////////////////////////
            """
            if ((i>0) and (Grid_states[i-1,j] == 0)) or ((i<length_x) and (Grid_states[i+1,j] == 0)):
    
                """
                Up and down
                """
                # Grid_states < 4 means there is Mo (3) or a Mo adatom (2)
                # We don't check S (1) or empty spaces (0)
                if (i>1) and (Grid_states[i-2,j] < 4): # Down
                    join_cluster_ij.append((i-2,j))    
                    join_sites += 1
                        
                if (i<length_x-1) and (Grid_states[i+2,j] < 4): # Up
                    join_cluster_ij.append((i+2,j))    
                    join_sites += 1
                        
                """
                Left up and down
                """
                if (j>1) and (i<length_x) and (Grid_states[i+1,j-2] < 4): # Left up
                    join_cluster_ij.append((i+1,j-2))    
                    join_sites += 1
                            
                if (j>1) and (i>0) and (Grid_states[i-1,j-2] < 4): # Left down
                    join_cluster_ij.append((i-1,j-2))    
                    join_sites += 1
                        
                """
                Right up and down
                """
                if (i<length_x) and (j<length_y) and (Grid_states[i+1,j+1] < 4): # Right up
                    join_cluster_ij.append((i+1,j+1))    
                    join_sites += 1
                                
                if (i>0) and (j<length_y) and (Grid_states[i-1,j+1] < 4): # Right down
                    join_cluster_ij.append((i-1,j+1))    
                    join_sites += 1
                        
                """
                //////////////////////////////// Right Mo ////////////////////////
                """
                
            if ((i>0) and (Grid_states[i-1,j] == 1)) or ((i<length_x) and (Grid_states[i+1,j] == 1)):
    
                """
                Up and down
                """
                if (i>1) and (Grid_states[i-2,j] < 4): # Down
                    join_cluster_ij.append((i-2,j))    
                    join_sites += 1
                        
                if (i<length_x-1) and (Grid_states[i+2,j] < 4): # Up
                    join_cluster_ij.append((i+2,j))    
                    join_sites += 1
                    
                """
                Left up and down
                """
                if (j>0) and (i<length_x) and (Grid_states[i+1,j-1] < 4): # Left up
                    join_cluster_ij.append((i+1,j-1))    
                    join_sites += 1
                        
                if (j>0) and (i>0) and (Grid_states[i-1,j-1] < 4): # Left down
                    join_cluster_ij.append((i-1,j-1))    
                    join_sites += 1
                        
                """
                Right up and down
                """
                if (i<length_x) and (j<length_y-1) and (Grid_states[i+1,j+2] < 4): # Right up
                    join_cluster_ij.append((i+1,j+2))    
                    join_sites += 1
                        
                if (i>0) and (j<length_y-1) and (Grid_states[i-1,j+2] < 4): # Right down
                    join_cluster_ij.append((i-1,j+2))    
                    join_sites += 1
                    
        
            if join_sites == 0: Grid_states[i,j] = 5 # Inner point of the crystal  
        # Grid points adjacent to the crystal, so they may join the crystal
        # Select unique elements from the list, sort them and transform into a list
        self.join_cluster_ij = sorted(join_cluster_ij) 
        self.Grid_states = Grid_states
                        
    def crystal_update(self,Grid_states,new_defect_ij,s,mig_defect):
        
        if (s+1 == 7) or (s+1 == 8):
            self.cluster_ij.append(new_defect_ij) # New defect incorporate to the crystal
            
            self.clustering_region(Grid_states,[new_defect_ij]) # Check new joining sites
            #join_cluster_ij = self.join_cluster_ij   
    
            # Remove elements of join_cluster_ij that already belongs to the cluster (cluster_ij) 
            self.join_cluster_ij = [x for x in self.cluster_ij+self.join_cluster_ij if x not in self.cluster_ij]
       
        elif (s > 9): # Migration of atoms at the crystal edge
            self.cluster_ij.remove(new_defect_ij) # This atoms changed his position
            self.cluster_ij.append(mig_defect) # New position of the atom
            
            self.join_cluster_ij = [] # There are no valid elements. We recalculate the region
            self.clustering_region(Grid_states,self.cluster_ij) # Check new joining sites
            # Remove elements of join_cluster_ij that already belongs to the cluster (cluster_ij) 
            self.join_cluster_ij = [x for x in self.cluster_ij+self.join_cluster_ij if x not in self.cluster_ij]
   
        # We update the Grid_states with inner points of the cluster (Grid_states = 5)
        return self.Grid_states
        