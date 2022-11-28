# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 20:15:19 2022

@author: ALDANADS
"""

from hex_lattice import*
import numpy as np
import matplotlib.pyplot as plt
from defects import Cluster
import random
from datetime import datetime
import shutil
import os 



def initialization(parameters,n_sim,save_data):

    # Random seed as time
    #random.seed(datetime.now())
    rng = np.random.default_rng() # Random Number Generator (RNG) object

    # Default resolution for figures
    plt.rcParams["figure.dpi"] = 300 # Default value of dpi = 300
    
    if save_data:
        files_copy = ['defects.py', 'hex_lattice.py', 'initialization.py','KMC.py','main_simulator.py']
        dst = r'C:\Users\aldanads\OneDrive - TCDUD.onmicrosoft.com\2D device simulator project\Publications\Layer growth\Simulations\Triangles_new\\'
        dst_data = save_simulation(files_copy,dst,n_sim) # Create folders and python files
    else:
        dst_data = ''
        
    
    """ --------------------------------------------------------------------------
     -----------------------Creating of the hexagonal grid -----------------------
     -----------------------------------------------------------------------------
    """ 
    a=0.639 # nm
    b=0.639 # nm
    device_size = (50, 50) # Size of the grid in nm. grid_size[0] is x and grid_size[1] is y.
    atom_colors=['orange','purple','blue', 'black'] # MoS2 -> First is Sulfur, second is Mo and third Vs
    
    # Activation energies --> Adatoms
    E_mig_zigzag = 1.2 # Zigzag direction
    E_mig_armchair = 1.2 # Armchair direction

    #E_nucleation = 1.7 # Kink nucleation (1.7 eV) --> Growing in armchair direction
    #E_propagation = 1.5 # Kink propagation (1.4 eV) --> Growing in zigzag direction
    E_nucleation = parameters[0]
    E_propagation = parameters[1]
    #E_desorption = 1.52
    E_desorption = 1.3
    # Activation energies --> Atoms at the crystal edge
    E_mig_armchair_edge = 3 # Armchair direction
    E_mig_zigzag_edge = 3 # Zigzag direction
    
    Backup_energy = [E_mig_zigzag,E_mig_armchair,E_nucleation,E_propagation,E_desorption,E_mig_zigzag_edge,E_mig_armchair_edge]
    
    # Temperature
    T = 1273
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
    MoS2_lattice = Hexagonal_lattice(a,b,device_size,atom_colors,Backup_energy,T)
    xv=MoS2_lattice.xv
    """
    distribution:
        'uniform' --> Uniform rows
        'skewed_gaussian' --> Skewed Gaussian distribution of defects
        'triangle' --> Right triangle distribution of defects
        'test 1: single adatom' --> Single adatom in the middle of the grid
    """
    distribution = ['uniform','skewed_gaussian','triangle','test 1: single adatom','Crystal seed']
    # skewness parameter --> a=0 is the normal distribution
    skewness = 12 # Skewness of the skewed Gaussian distribution
    fissure_region = (round(len(xv[0])/2)+2,100) # [0] middle point and [1] half width (nm)
    # Atomic specie -> Whay kind of atoms are affected by defects
    atomic_specie = 3 # Sulfur = 1 // Molibdenum = 3
    # Type of defects in the lattice
    defect_specie = 2 # Adatom = 2 // Crystal edge = 4 // Inner point of crystal = 5
    
    #prob_defects = adsortion_rate[n_sim]
    #prob_defects = 0.001
    prob_defects = parameters[2][n_sim]
    crystal_orientation = True
    pair_atom_defect=(3,4)
    MoS2_lattice.defect_distributions(prob_defects,fissure_region,skewness,distribution[4],pair_atom_defect,rng)
    pair_atom_defect = (atomic_specie,defect_specie)
    #MoS2_lattice.defect_distributions(prob_defects,fissure_region,skewness,distribution[0],pair_atom_defect,rng)
    
    distribution_parameters = [distribution[0],skewness,fissure_region,pair_atom_defect,prob_defects]

    
    MoS2_crystal = Cluster(MoS2_lattice.Grid_states) # Crystal seed
    
    MoS2_lattice.plot_lattice(crystal_orientation,'',0,0,True) # Initial state of the grid
    
    return MoS2_lattice,MoS2_crystal,distribution_parameters, dst_data,rng

def save_simulation(files_copy,dst,n_sim):
    

    parent_dir = 'Sim_'+str(n_sim)+'\\'
    os.makedirs(dst+parent_dir) 
    dst = dst+parent_dir
    program_directory = 'Program\\'
    data_directoy = 'Crystal growth\\'

    os.makedirs(dst + program_directory)
    os.makedirs(dst + data_directoy)

    dst_program = dst + program_directory

    for files in files_copy:
        shutil.copyfile(files, dst_program+files)
        
    return dst + data_directoy