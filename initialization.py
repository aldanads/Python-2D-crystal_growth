# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 20:15:19 2022

@author: ALDANADS
"""

from hex_lattice import*
import numpy as np
import matplotlib.pyplot as plt

def initialization():

    plt.rcParams["figure.dpi"] = 300 # Default value of dpi = 300
    
    """ --------------------------------------------------------------------------
     -----------------------Creating of the hexagonal grid -----------------------
     -----------------------------------------------------------------------------
    """ 
    a=0.639 # nm
    b=0.639 # nm
    device_size = (10 , 5) # Size of the grid in nm. grid_size[0] is x and grid_size[1] is y.
    atom_colors=['orange','purple','blue'] # MoS2 -> First is Sulfur, second is Mo and third Vs
    
    E_mig_armchair = 1
    E_mig_zigzag = 1
    
    Act_E = [E_mig_armchair,E_mig_zigzag]
    
    
    
    """  
    ---------------------------------------------
    --------------- Grid _states ----------------
    ---------------------------------------------
    Sulfur = 1
    Sulfur vacancy = 2
    Molibdenum = 3
    Empty space = 0
    -----------------------------------------------------------------------------
    """
    
    # Create MoS2 crystal
    MoS2_lattice = Hexagonal_lattice(a,b,device_size,atom_colors,Act_E)
    MoS2_lattice.create_hex_grid()
    xv=MoS2_lattice.xv
    yv=MoS2_lattice.yv
    """
    distribution:
        'uniform' --> Uniform rows
        'skewed_gaussian' --> Skewed Gaussian distribution of defects
        'triangle' --> Right triangle distribution of defects
    """
    
    
    distribution = 'uniform'
    prob_defects = 0.5 # prob of generating defects --> Peak density
    fissure_region = (round(len(xv[0])/2)+1,4) # [0] middle point and [1] half width (nm)
    
    # skewness parameter --> a=0 is the normal distribution
    skewness = 12 # Skewness of the skewed Gaussian distribution
    
    crystal_orientation = False
    MoS2_lattice.adam_atom()
    #MoS2_lattice.defect_distributions(prob_defects,fissure_region,skewness,distribution)
    MoS2_lattice.plot_lattice(crystal_orientation)
    
    return MoS2_lattice