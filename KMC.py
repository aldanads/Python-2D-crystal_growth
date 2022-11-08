# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:07:32 2022

@author: ALDANADS
"""
import numpy as np


def KMC(MoS2_lattice):
    
    Grid_states = MoS2_lattice.Grid_states
    coord_xy_Vs = np.where(Grid_states == 2)

    #allowed_events = events_available(Grid_states,coord_xy_Vs[0][0],coord_xy_Vs[1][0])
    allowed_events = events_available(Grid_states,4,30)
    return MoS2_lattice, allowed_events

def events_available(Grid_states,i,j):
    
    length_x = len(Grid_states)-1
    length_y = len(Grid_states[0])-1
    
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
    # 0: forbidden
    # 1: allowed
    # 2: left Mo
    # 3: right Mo
    

    if ((i>0) and (Grid_states[i-1,j] == 0)) or ((i<length_x) and (Grid_states[i+1,j] == 0)):
    # 0
    # 2: Mo (Grid_states) Left Mo --------------------------------------------
    # 0
    
        allowed_events[0] = 2 # Left Mo
        """
        Up and down
        """
        if (i>1) and Grid_states[i-2,j] == 3: # Down
           allowed_events[1] = 1 
           
        if (i<length_x-1) and (Grid_states[i+2,j] == 3): # Up
            allowed_events[2] = 1
            
        """
        Left up and down
        """
        if (j>1) and (i<length_x) and (Grid_states[i+1,j-2] == 3): # Left up
            allowed_events[3] = 1

        if (j>1) and (i>0) and (Grid_states[i-1,j-2] == 3): # Left down
            allowed_events[4] = 1
            
        """
        Right up and down
        """
        if (i<length_x) and (j<length_y) and (Grid_states[i+1,j+1] == 3): # Right up
            allowed_events[5] = 1
            
        if (i>0) and (j<length_y) and (Grid_states[i-1,j+1] == 3): # Right down
            allowed_events[6] = 1

            
    
    if ((i>0) and (Grid_states[i-1,j] == 1)) or ((i<length_x) and (Grid_states[i+1,j] == 1)):

    # S
    # 3: Mo (Grid_states) Right Mo --------------------------------------------
    # S
        allowed_events[0] = 3 # Right Mo
        
        """
        Up and down
        """
        if (i>1) and Grid_states[i-2,j] == 3: # Down
           allowed_events[1] = 1 
           
        if (i<length_x-1) and (Grid_states[i+2,j] == 3): # Up
            allowed_events[2] = 1

        """
        Left up and down
        """
        if (j>0) and (i<length_x) and (Grid_states[i+1,j-1] == 3): # Left up
            allowed_events[3] = 1
        
        if (j>0) and (i>0) and (Grid_states[i-1,j-1] == 3): # Left down
            allowed_events[4] = 1
            
        """
        Right up and down
        """
        if (i<length_x) and (j<length_y-1) and (Grid_states[i+1,j+2] == 3): # Right up
            allowed_events[5] = 1
                
        if (i>0) and (j<length_y-1) and (Grid_states[i-1,j+2] == 3): # Right down
            allowed_events[6] = 1
    
    
    
    return allowed_events

def TR(Ea,T):
    
    kb = 8.6173324E-5 # Boltzmann constant
    nu0=7E13;  # nu0 (s^-1) bond vibration frequency
    
    return nu0*np.exp(-Ea/(kb*T));