# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 14:03:55 2022

@author: ALDANADS
"""
import numpy as np

class Defects():
    
    
    def __init__(self,i,j,Backup_energy,atomic_specie,T,Grid_states,join_cluster_ij,split_regions):
        
        self.i = i
        self.j = j
        self.Backup_energy = Backup_energy

        self.energy_sets(i,j,Backup_energy,split_regions)
        
        self.atomic_specie = atomic_specie
        self.defect_specie = Grid_states[i,j]
        
        # Calculate the transition rates
        self.TR(Grid_states,join_cluster_ij)
        self.neighbors(i, j, Grid_states)
        
        
    """
    ---------------- Calculate neighbors -------------------------
    """
    
    def energy_sets(self,i,j,Backup_energy,split_regions):
        
        # Backup_energy[0] - zigzag
        # Backup_energy[1] - armchair
        # Backup_energy[2] - nucleation
        # Backup_energy[3] - propagation
        # Backup_energy[4] - desorption
        # Backup_energy[5] - E_mig_zigzag_edge
        # Backup_energy[6] - E_mig_armchair_edge
        
        if split_regions['key_parameter'] == 'adsorption rate':

            self.Act_E =[Backup_energy[0][0],Backup_energy[0][0],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1], # Migration at the substrate
                     Backup_energy[0][2], # Nucleation/Propagation
                     Backup_energy[0][4], # Desorption
                     Backup_energy[0][5],Backup_energy[0][5],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6]] # Migration at the edge
          
            self.TR_list= np.array([Backup_energy[1][0],Backup_energy[1][0],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1], # Migration substrate
                     Backup_energy[1][2], # Nucleation/Propagation
                     Backup_energy[1][4], # Desorption
                     Backup_energy[1][5],Backup_energy[1][5],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6]]) # Migration at the edge

        elif split_regions['key_parameter'] == 'migration substrate':
            
            if split_regions['Boundary'] == 'vertical right':
                
                if j > split_regions['Position']: 
                
                    self.Act_E =[Backup_energy[0][7],Backup_energy[0][7],Backup_energy[0][7],Backup_energy[0][7],Backup_energy[0][7],Backup_energy[0][7], # Migration at the substrate
                             Backup_energy[0][2], # Nucleation/Propagation
                             Backup_energy[0][4], # Desorption
                             Backup_energy[0][5],Backup_energy[0][5],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6]] # Migration at the edge
                  
                    self.TR_list= np.array([Backup_energy[1][7],Backup_energy[1][7],Backup_energy[1][7],Backup_energy[1][7],Backup_energy[1][7],Backup_energy[1][7], # Migration substrate
                             Backup_energy[1][2], # Nucleation/Propagation
                             Backup_energy[1][4], # Desorption
                             Backup_energy[1][5],Backup_energy[1][5],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6]]) # Migration at the edge
                else:
                    
                    self.Act_E =[Backup_energy[0][0],Backup_energy[0][0],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1], # Migration at the substrate
                             Backup_energy[0][2], # Nucleation/Propagation
                             Backup_energy[0][4], # Desorption
                             Backup_energy[0][5],Backup_energy[0][5],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6]] # Migration at the edge
                  
                    self.TR_list= np.array([Backup_energy[1][0],Backup_energy[1][0],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1], # Migration substrate
                             Backup_energy[1][2], # Nucleation/Propagation
                             Backup_energy[1][4], # Desorption
                             Backup_energy[1][5],Backup_energy[1][5],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6]]) # Migration at the edge

                
            elif split_regions['Boundary'] == 'vertical left':
                
                if j < split_regions['Position']: 
                    
                    self.Act_E =[Backup_energy[0][7],Backup_energy[0][7],Backup_energy[0][7],Backup_energy[0][7],Backup_energy[0][7],Backup_energy[0][7], # Migration at the substrate
                             Backup_energy[0][2], # Nucleation/Propagation
                             Backup_energy[0][4], # Desorption
                             Backup_energy[0][5],Backup_energy[0][5],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6]] # Migration at the edge
                 
                    self.TR_list= np.array([Backup_energy[1][7],Backup_energy[1][7],Backup_energy[1][7],Backup_energy[1][7],Backup_energy[1][7],Backup_energy[1][7], # Migration substrate
                             Backup_energy[1][2], # Nucleation/Propagation
                             Backup_energy[1][4], # Desorption
                             Backup_energy[1][5],Backup_energy[1][5],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6]]) # Migration at the edge
                else:
                    
                    self.Act_E =[Backup_energy[0][0],Backup_energy[0][0],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1], # Migration at the substrate
                             Backup_energy[0][2], # Nucleation/Propagation
                             Backup_energy[0][4], # Desorption
                             Backup_energy[0][5],Backup_energy[0][5],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6]] # Migration at the edge
      
                    self.TR_list= np.array([Backup_energy[1][0],Backup_energy[1][0],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1], # Migration substrate
                             Backup_energy[1][2], # Nucleation/Propagation
                             Backup_energy[1][4], # Desorption
                             Backup_energy[1][5],Backup_energy[1][5],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6]]) # Migration at the edge

            elif split_regions['Boundary'] == 'horizontal':
                
                if i > split_regions['Position']: 
                    self.Act_E =[Backup_energy[0][7],Backup_energy[0][7],Backup_energy[0][7],Backup_energy[0][7],Backup_energy[0][7],Backup_energy[0][7], # Migration at the substrate
                             Backup_energy[0][2], # Nucleation/Propagation
                             Backup_energy[0][4], # Desorption
                             Backup_energy[0][5],Backup_energy[0][5],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6]] # Migration at the edge
                  
                    self.TR_list= np.array([Backup_energy[1][7],Backup_energy[1][7],Backup_energy[1][7],Backup_energy[1][7],Backup_energy[1][7],Backup_energy[1][7], # Migration substrate
                             Backup_energy[1][2], # Nucleation/Propagation
                             Backup_energy[1][4], # Desorption
                             Backup_energy[1][5],Backup_energy[1][5],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6]]) # Migration at the edge
                else:
                    
                    self.Act_E =[Backup_energy[0][0],Backup_energy[0][0],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1], # Migration at the substrate
                             Backup_energy[0][2], # Nucleation/Propagation
                             Backup_energy[0][4], # Desorption
                             Backup_energy[0][5],Backup_energy[0][5],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6]] # Migration at the edge
      
                    self.TR_list= np.array([Backup_energy[1][0],Backup_energy[1][0],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1], # Migration substrate
                             Backup_energy[1][2], # Nucleation/Propagation
                             Backup_energy[1][4], # Desorption
                             Backup_energy[1][5],Backup_energy[1][5],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6]]) # Migration at the edge

        elif split_regions['key_parameter'] == 'desorption':
            
            if split_regions['Boundary'] == 'vertical right':
                
                if j > split_regions['Position']: 
                    self.Act_E =[Backup_energy[0][0],Backup_energy[0][0],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1], # Migration at the substrate
                             Backup_energy[0][2], # Nucleation/Propagation
                             Backup_energy[0][8], # Desorption
                             Backup_energy[0][5],Backup_energy[0][5],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6]] # Migration at the edge
                  
                    self.TR_list= np.array([Backup_energy[1][0],Backup_energy[1][0],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1], # Migration substrate
                             Backup_energy[1][2], # Nucleation/Propagation
                             Backup_energy[1][8], # Desorption
                             Backup_energy[1][5],Backup_energy[1][5],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6]]) # Migration at the edge

                else:
                    
                    self.Act_E =[Backup_energy[0][0],Backup_energy[0][0],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1], # Migration at the substrate
                             Backup_energy[0][2], # Nucleation/Propagation
                             Backup_energy[0][4], # Desorption
                             Backup_energy[0][5],Backup_energy[0][5],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6]] # Migration at the edge
                  
                    self.TR_list= np.array([Backup_energy[1][0],Backup_energy[1][0],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1], # Migration substrate
                             Backup_energy[1][2], # Nucleation/Propagation
                             Backup_energy[1][4], # Desorption
                             Backup_energy[1][5],Backup_energy[1][5],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6]]) # Migration at the edge

            if split_regions['Boundary'] == 'vertical left':
                
                if j < split_regions['Position']: 
                    self.Act_E =[Backup_energy[0][0],Backup_energy[0][0],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1], # Migration at the substrate
                             Backup_energy[0][2], # Nucleation/Propagation
                             Backup_energy[0][8], # Desorption
                             Backup_energy[0][5],Backup_energy[0][5],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6]] # Migration at the edge
                  
                    self.TR_list= np.array([Backup_energy[1][0],Backup_energy[1][0],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1], # Migration substrate
                             Backup_energy[1][2], # Nucleation/Propagation
                             Backup_energy[1][8], # Desorption
                             Backup_energy[1][5],Backup_energy[1][5],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6]]) # Migration at the edge

                else:
                    
                    self.Act_E =[Backup_energy[0][0],Backup_energy[0][0],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1], # Migration at the substrate
                             Backup_energy[0][2], # Nucleation/Propagation
                             Backup_energy[0][4], # Desorption
                             Backup_energy[0][5],Backup_energy[0][5],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6]] # Migration at the edge
                  
                    self.TR_list= np.array([Backup_energy[1][0],Backup_energy[1][0],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1], # Migration substrate
                             Backup_energy[1][2], # Nucleation/Propagation
                             Backup_energy[1][4], # Desorption
                             Backup_energy[1][5],Backup_energy[1][5],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6]]) # Migration at the edge

            if split_regions['Boundary'] == 'horizontal':
                
                if i > split_regions['Position']:  
                    self.Act_E =[Backup_energy[0][0],Backup_energy[0][0],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1], # Migration at the substrate
                             Backup_energy[0][2], # Nucleation/Propagation
                             Backup_energy[0][8], # Desorption
                             Backup_energy[0][5],Backup_energy[0][5],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6]] # Migration at the edge
                  
                    self.TR_list= np.array([Backup_energy[1][0],Backup_energy[1][0],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1], # Migration substrate
                             Backup_energy[1][2], # Nucleation/Propagation
                             Backup_energy[1][8], # Desorption
                             Backup_energy[1][5],Backup_energy[1][5],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6]]) # Migration at the edge

                else:
                    
                    self.Act_E =[Backup_energy[0][0],Backup_energy[0][0],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1],Backup_energy[0][1], # Migration at the substrate
                             Backup_energy[0][2], # Nucleation/Propagation
                             Backup_energy[0][4], # Desorption
                             Backup_energy[0][5],Backup_energy[0][5],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6],Backup_energy[0][6]] # Migration at the edge
                  
                    self.TR_list= np.array([Backup_energy[1][0],Backup_energy[1][0],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1],Backup_energy[1][1], # Migration substrate
                             Backup_energy[1][2], # Nucleation/Propagation
                             Backup_energy[1][4], # Desorption
                             Backup_energy[1][5],Backup_energy[1][5],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6],Backup_energy[1][6]]) # Migration at the edge

    
    
    def neighbors(self,i,j,Grid_states):
        
        List_neighbors = []
        length_x = len(Grid_states)-1
        length_y = len(Grid_states[0])-1
        
        if ((i>0) and (Grid_states[i-1,j] == 0)) or ((i<length_x) and (Grid_states[i+1,j] == 0)):
        # 0
        # 2: Mo (Grid_states) Left Mo --------------------------------------------
        # 0
        
            """
            Check up and down
            """
            if (i>1) and ((Grid_states[i-2,j] == 2) or (Grid_states[i-2,j] == 4)): # Down
               List_neighbors.append((i-2,j))
               
            if (i<length_x-1) and ((Grid_states[i+2,j] == 2) or (Grid_states[i+2,j] == 4)): # Up
               List_neighbors.append((i+2,j))
               
            """
            Left up and down
            """
            if (j>1) and (i<length_x) and ((Grid_states[i+1,j-2] == 2) or (Grid_states[i+1,j-2] == 4)): # Left up
                List_neighbors.append((i+1,j-2))

            if (j>1) and (i>0) and ((Grid_states[i-1,j-2] == 2) or (Grid_states[i-1,j-2] == 4)): # Left down
                List_neighbors.append((i-1,j-2))
                
            """
            Right up and down
            """
            if (i<length_x) and (j<length_y) and ((Grid_states[i+1,j+1] == 2) or (Grid_states[i+1,j+1] == 4)): # Right up
                List_neighbors.append((i+1,j+1))
                   
            if (i>0) and (j<length_y) and ((Grid_states[i-1,j+1] == 2) or (Grid_states[i-1,j+1] ==  4)): # Right down
                List_neighbors.append((i-1,j+1))
                
        if ((i>0) and (Grid_states[i-1,j] == 1)) or ((i<length_x) and (Grid_states[i+1,j] == 1)):
        # S
        # 3: Mo (Grid_states) Right Mo --------------------------------------------
        # S
 
            """
            Up and down
            """
            if (i>1) and (Grid_states[i-2,j] == 2 or (Grid_states[i-2,j] == 4)): # Down
                List_neighbors.append((i-2,j)) 
               
            if (i<length_x-1) and ((Grid_states[i+2,j] == 2) or (Grid_states[i+2,j] == 4)): # Up
                List_neighbors.append((i+2,j)) 

            """
            Left up and down
            """
            if (j>0) and (i<length_x) and ((Grid_states[i+1,j-1] == 2) or (Grid_states[i+1,j-1] == 4)): # Left up
                List_neighbors.append((i+1,j-1)) 
            
            if (j>0) and (i>0) and ((Grid_states[i-1,j-1] == 2) or (Grid_states[i-1,j-1] == 4)): # Left down
                List_neighbors.append((i-1,j-1)) 
                
            """
            Right up and down
            """
            if (i<length_x) and (j<length_y-1) and ((Grid_states[i+1,j+2] == 2) or (Grid_states[i+1,j+2] == 4)): # Right up
                List_neighbors.append((i+1,j+2)) 

                    
            if (i>0) and (j<length_y-1) and ((Grid_states[i-1,j+2] == 2) or (Grid_states[i-1,j+2] == 4)): # Right down
                List_neighbors.append((i-1,j+2)) 

        self.List_neighbors = List_neighbors
        return List_neighbors

    
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
            # Desorption - allowed_events[8]
            # Migration edge - allowed_events[9-14]
        """

        self.allowed_events = np.zeros(len(self.Act_E)+1)

        # Adatoms --> Migrate, nucleate or desorpt
        if Grid_states[self.i,self.j] == 2: 
    
            #self.allowed_events = allowed_events 
            # allowed_events[0] =
            # 2: left Mo
            # 3: right Mo
            # allowed_events[1:end] =
            # 0: forbidden
            # 1: allowed
            
            self.mig_available(Grid_states) # Migration

            self.nucleation_available(join_cluster_ij,Grid_states) # Nucleation/propagation
            
            self.allowed_events[8] = 1
            
        # Atoms at the edge of the crystal --> They can migrate through the edge
        if Grid_states[self.i,self.j] == 4:
            
            #self.allowed_events = allowed_events 
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
            
    # Calculate the nucleation region
    def nucleation_available(self,join_cluster_ij,Grid_states):
        
        i = self.i
        j = self.j
        length_x = len(Grid_states)-1
        length_y = len(Grid_states[0])-1

        join_cluster_ij = list(set(join_cluster_ij))
 
        # Defect in (i,j) is in contact with the cluster
        if (i,j) in join_cluster_ij: 
            
            # Kink propagation --> Zigzag
            if (((i>1) and Grid_states[i-2,j]) >= 4) or ((i<length_x-1) and (Grid_states[i+2,j] >= 4)):
                self.allowed_events[7] = 1
                self.Act_E[6] = self.Backup_energy[0][3] # Set the nucleation energy as zigzag
                return  # If it is part of one zigzag, it is enough
                
            # Kink nucleation --> Armchair
            # Search in a square 3x4 around (i,j) --> Remember that init:final:steps exclude final
            #if (i>0 and i < length_x and j > 1 and j < length_y-2) and (Grid_states[i-1:i+2:1,j-2:j+3:1] == 4).any():
            #    self.allowed_events[7] = 1
            
            
            
            
            # Kink nucleation --> Armchair
            # We check if the position (i,j) have any crystal in contact
            # Then we have to check if the position (i,j) follow a zigzag direction in diagonal
            if ((i>0) and (Grid_states[i-1,j] == 0)) or ((i<length_x) and (Grid_states[i+1,j] == 0)):
            # 0
            # 2: Mo (Grid_states) Left Mo --------------------------------------------
            # 0 
                """
                Check if there is a crystal at left up and down
                """
                
                if (j > 1) and (i < length_x) and (Grid_states[i+1,j-2] >= 4): # Left up
                    self.allowed_events[7] = 1
                    if (j > 2) and (i < length_x - 1) and Grid_states[i+2,j-3] >= 4: # The continuation of a diagonal zigzag row
                        self.Act_E[6] = self.Backup_energy[0][3] # Set the nucleation energy as zigzag
                        return # If it is part of one zigzag row, it is enough
                    
    
                if (j > 1) and (i > 0) and (Grid_states[i-1,j-2] >= 4): # Left down
                    self.allowed_events[7] = 1
                    if (j > 2) and (i > 1) and Grid_states[i-2,j-3] >= 4: # The continuation of a diagonal zigzag row
                        self.Act_E[6] = self.Backup_energy[0][3] # Set the nucleation energy as zigzag
                        return # If it is part of one zigzag row, it is enough
                
                    
                """
                Check if there is a crystal at right up and down
                """
                
                if (i < length_x) and (j < length_y) and (Grid_states[i+1,j+1] >= 4): # Right up
                    self.allowed_events[7] = 1
                    if (i < length_x - 1) and (j < length_y-2) and Grid_states[i+2,j+3] >= 4: # The continuation of a diagonal zigzag row
                        self.Act_E[6] = self.Backup_energy[0][3] # Set the nucleation energy as zigzag
                        return # If it is part of one zigzag row, it is enough
                    
                if (i > 0) and (j < length_y) and (Grid_states[i-1,j+1] >= 4): # Right down
                    self.allowed_events[7] = 1
                    if (i > 1) and (j < length_y-2) and Grid_states[i-2,j+3] >= 4: # The continuation of a diagonal zigzag row
                        self.Act_E[6] = self.Backup_energy[0][3] # Set the nucleation energy as zigzag
                        return # If it is part of one zigzag row, it is enough
                    
            
            
            if ((i>0) and (Grid_states[i-1,j] == 1)) or ((i<length_x) and (Grid_states[i+1,j] == 1)):
            # S
            # 3: Mo (Grid_states) Right Mo --------------------------------------------
            # S
            
            
                """
                Check if there is a crystal at left up and down
                """
                
                if (j>0) and (i < length_x) and (Grid_states[i+1,j-1] >= 4): # Left up
                    self.allowed_events[7] = 1
                    if (j>2) and (i < length_x - 1) and Grid_states[i+2,j-3] >= 4: # The continuation of a diagonal zigzag row
                        self.Act_E[6] = self.Backup_energy[0][3] # Set the nucleation energy as zigzag
                        return # If it is part of one zigzag row, it is enough
                
                if (j>0) and (i>0) and (Grid_states[i-1,j-1] >= 4): # Left down
                    self.allowed_events[7] = 1
                    if (j > 2) and (i > 1) and Grid_states[i-2,j-3] >= 4: # The continuation of a diagonal zigzag row
                        self.Act_E[6] = self.Backup_energy[0][3] # Set the nucleation energy as zigzag
                        return # If it is part of one zigzag row, it is enough
                    
                
                    
                """
                Check if there is a crystal at right up and down
                """
                
                if (i<length_x) and (j<length_y-1) and (Grid_states[i+1,j+2] >= 4): # Right up
                    self.allowed_events[7] = 1
                    if (i < length_x - 1) and (j < length_y-2) and Grid_states[i+2,j+3] >= 4: # The continuation of a diagonal zigzag row
                        self.Act_E[6] = self.Backup_energy[0][3] # Set the nucleation energy as zigzag
                        return # If it is part of one zigzag row, it is enough
                        
                if (i > 0) and (j < length_y-1) and (Grid_states[i-1,j+2] >= 4): # Right down
                    self.allowed_events[7] = 1
                    if (i > 1) and (j < length_y-2) and Grid_states[i-2,j+3] >= 4: # The continuation of a diagonal zigzag row
                        self.Act_E[6] = self.Backup_energy[0][3] # Set the nucleation energy as zigzag
                        return # If it is part of one zigzag row, it is enough
                    
                
            
                
            
    # Migration of atoms at the edge of the crystal
    def migration_edge(self,join_cluster_ij,Grid_states):

        i = self.i
        j = self.j
        
        # If the atom is surrounded by more than 2 other atoms, it is immobile
        # Minus 1 because we don't take into account the defect at the center (i,j)
        #if type(Grid_states) != np.ndarray:
        #if (i>1) and sum(sum(Grid_states[i-2:i+3:1,j-2:j+3:1] >= 4))-1 > 2: return
        if i<2: return
        
        # In case the particle form a chain --> Have particles 4 and 5 in opposite sites 
        if (4 in Grid_states[i-2:i:1,j-2:j+3:1] and 5 in Grid_states[i+1:i+3:1,j-2:j+3:1]):
            return
        elif (5 in Grid_states[i-2:i:1,j-2:j+3:1] and 4 in Grid_states[i+1:i+3:1,j-2:j+3:1]): 
            return
        elif (4 in Grid_states[i-2:i+3:1,j-2:j:1] and 5 in Grid_states[i-2:i+3:1,j+1:j+3:1]):
            return
        elif (5 in Grid_states[i-2:i+3:1,j-2:j:1] and 4 in Grid_states[i-2:i+3:1,j+1:j+3:1]):
            return
        length_x = len(Grid_states)-1
        length_y = len(Grid_states[0])-1
        atomic_specie = self.atomic_specie

    
        # We select the position supported by at least two atoms of the crystal
        # That is, the defect in the edge can't jump outside the crystal during the migration process
        # Jump outside the crystal is the detach process
        join_cluster_ij = [x for x in join_cluster_ij if (join_cluster_ij.count(x) > 2) or (join_cluster_ij.count(x) == 2 and sum(sum(Grid_states[i-2:i+3:1,j-2:j+3:1] >= 4))-1 >= 1)]
        #                    or (join_cluster_ij.count(x) == 2 and sum(sum(Grid_states[i-2:i+3:1,j-2:j+3:1] >= 4))-1 == ]
        
        """
        Up and down
        """
        if (i>1) and (((i-2,j) in join_cluster_ij) and Grid_states[i-2,j] == atomic_specie): # Down
            self.allowed_events[9] = 1 
                   
        if (i<length_x-1) and (((i+2,j) in join_cluster_ij) and Grid_states[i+2,j] == atomic_specie): # Up
            self.allowed_events[10] = 1
            
            
        
        if ((i>0) and (Grid_states[i-1,j] == 0)) or ((i<length_x) and (Grid_states[i+1,j] == 0)):
        # 0
        # 2: Mo (Grid_states) Left Mo --------------------------------------------
        # 0
                
            self.allowed_events[0] = 2 # Left Mo

                
            """
            Left up and down
            """
            if (j>1) and (i<length_x) and (((i+1,j-2) in join_cluster_ij) and Grid_states[i+1,j-2] == atomic_specie): # Left up
                self.allowed_events[11] = 1

            if (j>1) and (i>0) and (((i-1,j-2) in join_cluster_ij) and Grid_states[i-1,j-2] == atomic_specie): # Left down
                self.allowed_events[12] = 1
                
            """
            Right up and down
            """
            if (i<length_x) and (j<length_y) and (((i+1,j+1) in join_cluster_ij) and Grid_states[i+1,j+1] == atomic_specie): # Right up
                self.allowed_events[13] = 1
                
            if (i>0) and (j<length_y) and (((i-1,j+1) in join_cluster_ij) and Grid_states[i-1,j+1] == atomic_specie): # Right down
                self.allowed_events[14] = 1
                
                

        if ((i>0) and (Grid_states[i-1,j] == 1)) or ((i<length_x) and (Grid_states[i+1,j] == 1)):
        # S
        # 3: Mo (Grid_states) Right Mo --------------------------------------------
        # S
        
            self.allowed_events[0] = 3 # Right Mo
            
            
            """
            Left up and down
            """
            if (j>0) and (i<length_x) and (((i+1,j-1) in join_cluster_ij) and Grid_states[i+1,j-1] == atomic_specie): # Left up
                self.allowed_events[11] = 1
            
            if (j>0) and (i>0) and (((i-1,j-1) in join_cluster_ij) and Grid_states[i-1,j-1] == atomic_specie): # Left down
                self.allowed_events[12] = 1
                
            """
            Right up and down
            """
            if (i<length_x) and (j<length_y-1) and (((i+1,j+2) in join_cluster_ij) and Grid_states[i+1,j+2] == atomic_specie): # Right up
                self.allowed_events[13] = 1
                    
            if (i>0) and (j<length_y-1) and (((i-1,j+2) in join_cluster_ij) and Grid_states[i-1,j+2] == atomic_specie): # Right down
                self.allowed_events[14] = 1
        
        
        """
        ________________________ TRANSITION RATES __________________________
        """
        
        
    def TR(self,Grid_states,join_cluster_ij):
        
        self.events_available(Grid_states,join_cluster_ij)
        #kb = 8.6173324E-5 # Boltzmann constant
        #nu0=7E13;  # nu0 (s^-1) bond vibration frequency
        #nu0 = 1E12 #  Shuai, Chen, Gao Junfeng, Bharathi M. Srinivasan, and Zhang Yong-Wei. "A kinetic Monte Carlo study for mono-and bi-layer growth of MoS2 during chemical vapor deposition." Acta Physico-Chimica Sinica 35, no. 10 (2019): 1119-1127.
        allowed_events = self.allowed_events
        TR = np.zeros(len(allowed_events)-1)
           
        #TR = nu0*np.exp(-np.array(self.Act_E)/(kb*T))
        print(len(TR))
        print(len(self.TR_list))
        print(len(allowed_events))
        TR[allowed_events[1:] != 0] = self.TR_list[allowed_events[1:] != 0]
        
        if self.Act_E[6] == self.Backup_energy[0][3]:
            TR[6] = self.Backup_energy[1][3]
        self.TR = TR
        
        return TR
        
        """
        ________________________ PROCESSES __________________________
        """
        

    def processes(self,Grid_states,s,i,j):
        
        s = s+1 # s is selected from TR, which is smaller than allowed_events
        
        # i = self.i
        # j = self.j
        
        atomic_specie = self.atomic_specie
        defect_specie = self.defect_specie
        length_x = len(Grid_states)-1

        
        """ Adatom migration: Mo migrating to Mo positions
            # Mo left (2) or Mo right (3) allowed_events[0]
            # Down - Zigzag - allowed_events[1]    
            # Up - Zigzag - allowed_events[2]   
            # Left up - Armchair - allowed_events[3]
            # Left down - Armchair - allowed_events[4]
            # Right up - Armchair - allowed_events[5]
            # Right down - Armchair - allowed_events[6]
            # Nucleation - allowed_events[7]
            # Desorption - allowed_events[8]
            # Migration edge - allowed_events[9-14]
        """
        if (s == 7): # The defect join the cluster -> Exit the function
            Grid_states[i,j] = 4 
            return Grid_states
        
        if (s == 8): # Desorption --> Exit the function 
            Grid_states[i,j] = 3
            return Grid_states
        
        # Mo right (2)
        elif ((i>0) and (Grid_states[i-1,j] == 0)) or ((i<length_x) and (Grid_states[i+1,j] == 0)): 
        
            """
            ---------------------- Migration ----------------------------------
            """
        
            if (s == 1 or s == 9): # Down - Zigzag
                Grid_states[i,j] = atomic_specie
                Grid_states[i-2,j] = defect_specie
                self.i = i-2
                
            elif (s == 2 or s == 10): # Up - Zigzag
                Grid_states[i,j] = atomic_specie
                Grid_states[i+2,j] = defect_specie
                self.i = i+2
                
            elif (s == 3 or s == 11): # Left up - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i+1,j-2] = defect_specie
                self.i = i+1
                self.j = j-2
            
            elif (s == 4 or s == 12): # Left down - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i-1,j-2] = defect_specie
                self.i = i-1
                self.j = j-2
                
            elif (s == 5 or s == 13): # Right up - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i+1,j+1] = defect_specie
                self.i = i+1
                self.j = j+1
                
            elif (s == 6 or s == 14): # Right down - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i-1,j+1] = defect_specie
                self.i = i-1
                self.j = j+1
        
        # Mo right (3)
        elif ((i>0) and (Grid_states[i-1,j] == 1)) or ((i<length_x) and (Grid_states[i+1,j] == 1)): 
        
            """
            ---------------------- Migration ----------------------------------
            """

            
            if (s == 1 or s == 9): # Down - Zigzag
                Grid_states[i,j] = atomic_specie
                Grid_states[i-2,j] = defect_specie
                self.i = i-2
                
            elif (s == 2 or s == 10): # Up - Zigzag
                Grid_states[i,j] = atomic_specie
                Grid_states[i+2,j] = defect_specie
                self.i = i+2
                
            elif (s == 3 or s == 11): # Left up - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i+1,j-1] = defect_specie
                self.i = i+1
                self.j = j-1
            
            elif (s == 4 or s == 12): # Left down - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i-1,j-1] = defect_specie
                self.i = i-1
                self.j = j-1
            
            elif (s == 5 or s == 13): # Right up - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i+1,j+2] = defect_specie
                self.i = i+1
                self.j = j+2
            
            elif (s == 6 or s == 14): # Right down - Armchair
                Grid_states[i,j] = atomic_specie
                Grid_states[i-1,j+2] = defect_specie
                self.i = i-1
                self.j = j+2
            
        return Grid_states
    




            
class Cluster():
        
    def __init__(self,Grid_states):
            
        #Calculate the crystal coordinates
        cluster_ij = np.where(Grid_states >= 4) # Coordinates of the full crystal
        cluster_ij = [(cluster_ij[0][i],cluster_ij[1][i]) for i in np.arange(len(cluster_ij[0]))]
        self.cluster_ij = cluster_ij # List of tuples with the coordinates of the crystal

        self.cluster_size = sum(sum(Grid_states == 4)) + sum(sum(Grid_states == 5))
        self.v_cluster_size = [self.cluster_size]
        self.Mo_sites = sum(sum(Grid_states > 1))
        self.coverage = self.cluster_size / self.Mo_sites
        
        # Create a empty list --> Add tuples with coordinates where defects can 
        # joining to the cluster
        self.join_cluster_ij = []   
        
        # Calculate the region where adatoms can join the crystal
        # Update: Inner part of the crystal is 5
        # The edge remains 4
        self.clustering_region(Grid_states,cluster_ij) 
        

        edge_ij = np.where(self.Grid_states == 4) # Coordinates of the edge
        edge_ij = [(edge_ij[0][i],edge_ij[1][i]) for i in np.arange(len(edge_ij[0]))] # Convert to tuples
        self.edge_ij = edge_ij # List of tuples with the coordinates of the edge


            
        
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
            triangle = np.zeros(3, dtype=int)+2 # Borders of the domain turn Grid_states[i,j] = 5

                
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
                    triangle[0] -= 1
                        
                if (i<length_x-1) and (Grid_states[i+2,j] < 4): # Up
                    join_cluster_ij.append((i+2,j))    
                    join_sites += 1
                    triangle[1] -= 1

                        
                """
                Left up and down
                """
                if (j>1) and (i<length_x) and (Grid_states[i+1,j-2] < 4): # Left up
                    join_cluster_ij.append((i+1,j-2))    
                    join_sites += 1
                    triangle[2] -= 1

                            
                if (j>1) and (i>0) and (Grid_states[i-1,j-2] < 4): # Left down
                    join_cluster_ij.append((i-1,j-2))    
                    join_sites += 1
                    triangle[2] -= 1

                        
                """
                Right up and down
                """
                if (i<length_x) and (j<length_y) and (Grid_states[i+1,j+1] < 4): # Right up
                    join_cluster_ij.append((i+1,j+1))    
                    join_sites += 1
                    triangle[1] -= 1

                                
                if (i>0) and (j<length_y) and (Grid_states[i-1,j+1] < 4): # Right down
                    join_cluster_ij.append((i-1,j+1))    
                    join_sites += 1
                    triangle[0] -= 1

                        
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
                    triangle[0] -= 1

                        
                if (i<length_x-1) and (Grid_states[i+2,j] < 4): # Up
                    join_cluster_ij.append((i+2,j))    
                    join_sites += 1
                    triangle[1] -= 1
                    
                """
                Left up and down
                """
                if (j>0) and (i<length_x) and (Grid_states[i+1,j-1] < 4): # Left up
                    join_cluster_ij.append((i+1,j-1))    
                    join_sites += 1
                    triangle[2] -= 1
                        
                if (j>0) and (i>0) and (Grid_states[i-1,j-1] < 4): # Left down
                    join_cluster_ij.append((i-1,j-1))    
                    join_sites += 1
                    triangle[2] -= 1
                        
                """
                Right up and down
                """
                if (i<length_x) and (j<length_y-1) and (Grid_states[i+1,j+2] < 4): # Right up
                    join_cluster_ij.append((i+1,j+2))    
                    join_sites += 1
                    triangle[1] -= 1
                        
                if (i>0) and (j<length_y-1) and (Grid_states[i-1,j+2] < 4): # Right down
                    join_cluster_ij.append((i-1,j+2))    
                    join_sites += 1
                    triangle[0] -= 1
        
            #if join_sites == 0: Grid_states[i,j] = 5 # Inner point of the crystal
            if any(triangle == 2): Grid_states[i,j] = 5 # Form a triangle
        # Grid points adjacent to the crystal, so they may join the crystal
        # Select unique elements from the list, sort them and transform into a list
        self.join_cluster_ij = join_cluster_ij 
        self.Grid_states = Grid_states
                        
    def crystal_update(self,Grid_states,new_defect_ij,s,mig_defect):
        
        
        """ 
            # Nucleation - allowed_events[7]
            # Desorption - allowed_events[8]
            # Migration edge - allowed_events[9-14]
        """
        if (s+1 == 7): # Nucleation
            self.cluster_ij.append(new_defect_ij) # New defect incorporate to the crystal
            
            self.clustering_region(Grid_states,[new_defect_ij]) # Check new joining sites
            if self.Grid_states[mig_defect] == 4: self.edge_ij.append(mig_defect) #Update the edge
            self.cluster_size += 1
            self.coverage = self.cluster_size / self.Mo_sites # The coverage of the layer by the crystal          
            # Remove elements of join_cluster_ij that already belongs to the cluster (edge_ij) 
            self.join_cluster_ij = [x for x in self.cluster_ij+self.join_cluster_ij if x not in self.cluster_ij]
       

        elif (s+1 > 8): # Migration of atoms at the crystal edge
            self.cluster_ij.remove(new_defect_ij) # This atoms changed his position
            self.cluster_ij.append(mig_defect) # New position of the atom

            self.join_cluster_ij = [] # There are no valid elements. We recalculate the region
            self.clustering_region(Grid_states,self.cluster_ij) # Check new joining sites
            self.edge_ij.remove(new_defect_ij)
            if self.Grid_states[mig_defect] == 4: self.edge_ij.append(mig_defect) # Update the edge

            # Remove elements of join_cluster_ij that already belongs to the cluster (edge_ij) 
            self.join_cluster_ij = [x for x in self.cluster_ij+self.join_cluster_ij if x not in self.cluster_ij]
   
        # We update the Grid_states with inner points of the cluster (Grid_states = 5)
        return self.Grid_states
    
    def crystal_area(self):
        
        self.v_cluster_size.append(self.cluster_size)

        
       