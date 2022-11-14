# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 20:15:19 2022

@author: ALDANADS
"""

from hex_lattice import*
import numpy as np
import matplotlib.pyplot as plt
from defects import Cluster


def initialization():

    plt.rcParams["figure.dpi"] = 300 # Default value of dpi = 300
    
    """ --------------------------------------------------------------------------
     -----------------------Creating of the hexagonal grid -----------------------
     -----------------------------------------------------------------------------
    """ 
    a=0.639 # nm
    b=0.639 # nm
    device_size = (50, 50) # Size of the grid in nm. grid_size[0] is x and grid_size[1] is y.
    atom_colors=['orange','purple','blue', 'black'] # MoS2 -> First is Sulfur, second is Mo and third Vs
    
    # Activation energies
    E_mig_armchair = 1 # Armchair direction
    E_mig_zigzag = 1 # Zigzag direction
    E_nucleation = 1.025 # Kink nucleation --> Growing in armchair direction
    E_propagation = 0.85 # Kink propagation --> Growing in zigzag direction

    
    Act_E = [E_mig_zigzag,E_mig_zigzag,E_mig_armchair,E_mig_armchair,E_mig_armchair,E_mig_armchair,E_nucleation,E_propagation]
    
    # Temperature
    T = 300
    # T = 835 Celsius in exp
    # time for growing: 5 min
    
    """  
    ---------------------------------------------
    --------------- Grid _states ----------------
    ---------------------------------------------
    Sulfur = 1
    Defect = 2
    Molibdenum = 3
    Empty space = 0
    -----------------------------------------------------------------------------
    """
    
    # Create MoS2 crystal
    MoS2_lattice = Hexagonal_lattice(a,b,device_size,atom_colors,Act_E,T)
    xv=MoS2_lattice.xv
    """
    distribution:
        'uniform' --> Uniform rows
        'skewed_gaussian' --> Skewed Gaussian distribution of defects
        'triangle' --> Right triangle distribution of defects
        'test 1: single adatom' --> Single adatom in the middle of the grid
    """
    distribution = ['uniform','skewed_gaussian','triangle','test 1: single adatom']
    # skewness parameter --> a=0 is the normal distribution
    skewness = 12 # Skewness of the skewed Gaussian distribution
    fissure_region = (round(len(xv[0])/2)+8,6) # [0] middle point and [1] half width (nm)
    # Atomic specie -> Whay kind of atoms are affected by defects
    atomic_specie = 3 # Sulfur = 1 // Molibdenum = 3
    # Type of defects in the lattice
    defect_specie = 2 # Adatom = 2 // Crystal edge = 4 // Inner point of crystal = 5
    
    prob_defects = 0.25 # prob of generating defects --> Peak density
    

    crystal_orientation = False
    pair_atom_defect=(3,4)
    MoS2_lattice.defect_distributions(prob_defects,fissure_region,skewness,distribution[3],pair_atom_defect)
    pair_atom_defect = (atomic_specie,defect_specie)
    MoS2_lattice.defect_distributions(prob_defects,fissure_region,skewness,distribution[0],pair_atom_defect)

    
    MoS2_crystal = Cluster(MoS2_lattice.Grid_states) # Crystal seed
    
    MoS2_lattice.plot_lattice(crystal_orientation) # Initial state of the grid
    
    return MoS2_lattice,MoS2_crystal,pair_atom_defect