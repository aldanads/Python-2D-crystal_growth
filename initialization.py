# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 20:15:19 2022

@author: ALDANADS
"""

from hex_lattice import Hexagonal_lattice
import numpy as np
import matplotlib.pyplot as plt
from defects import Cluster
#import random
#from datetime import datetime
import shutil
import os 



def initialization(parameters,n_sim,save_data):

    # Random seed as time
    rng = np.random.default_rng() # Random Number Generator (RNG) object

    # Default resolution for figures
    plt.rcParams["figure.dpi"] = 300 # Default value of dpi = 300
    
    if save_data:
        files_copy = ['defects.py', 'hex_lattice.py', 'initialization.py','KMC.py','main_simulator.py','load_variables.py']
        dst = r'C:\Users\aldanads\OneDrive - TCDUD.onmicrosoft.com\2D device simulator project\Publications\Layer growth\Simulations\Triangles_new\Adsorption_rate_6\\'
        paths = save_simulation(files_copy,dst,n_sim) # Create folders and python files
    else:
        paths = {'data': ''}
        
    
    """ --------------------------------------------------------------------------
     -----------------------Creating of the hexagonal grid -----------------------
     -----------------------------------------------------------------------------
    """ 
    a=0.639 # nm
    b=0.639 # nm
    device_size = (50, 50) # Size of the grid in nm. grid_size[0] is x and grid_size[1] is y.
    atom_colors=['orange','purple','blue', 'black'] # MoS2 -> First is Sulfur, second is Mo and third Vs
    
# =============================================================================
#     # Activation energies 
# =============================================================================
    E_mig_zigzag = 1.3 # Zigzag direction
    E_mig_armchair = 1.3 # Armchair direction

    E_nucleation = 1.7 # Kink nucleation (1.7 eV) --> Growing in armchair direction
    E_propagation = 1.4 # Kink propagation (1.4 eV) --> Growing in zigzag direction
    #E_desorption = 1.52
    E_desorption = 1.4
 
    
# =============================================================================
#     # Edge diffusion comparable to the attachment barrier
#     # Edge diffusion = 1.84 eV
#     Li, Xiaying, Shiping Zhang, Shuai Chen, Xingli Zhang, Junfeng Gao, Yong-Wei Zhang, Jijun Zhao et al. 
#     "Mo concentration controls the morphological transitions from dendritic to semicompact, and to compact growth of monolayer crystalline MoS2 on various substrates." 
#     ACS applied materials & interfaces 11, no. 45 (2019): 42751-42759.
# =============================================================================

    E_mig_armchair_edge = 1.84 # Armchair direction
    E_mig_zigzag_edge = 1.84 # Zigzag direction
    
    Backup_energy = [E_mig_zigzag,E_mig_armchair,E_nucleation,E_propagation,E_desorption,E_mig_zigzag_edge,E_mig_armchair_edge]
    
    # Temperature
    T = 1108
    #T = parameters[2][n_sim]
    # T = 835 Celsius in exp
    # time for growing: 5 min
    
    """  
    ---------------------------------------------
    --------------- Grid _states ----------------
    ---------------------------------------------
    Sulfur = 1
    Mo adatom = 2
    Molibdenum = 3
    Atom at the edge of the crystal = 4
    Atom part of the crystal = 5
    Empty space = 0
    -----------------------------------------------------------------------------
    """
    
    # Create MoS2 crystal
    MoS2_lattice = Hexagonal_lattice(a,b,device_size,atom_colors,Backup_energy,T)
    xv=MoS2_lattice.xv
    """
# =============================================================================
#     distribution:
#         'uniform' --> Uniform rows
#         'skewed_gaussian' --> Skewed Gaussian distribution of defects
#         'triangle' --> Right triangle distribution of defects
#         'test 1: single adatom' --> Single adatom in the middle of the grid
# =============================================================================
    """
    distribution = ['uniform','skewed_gaussian','triangle','test 1: single adatom','test 2: column defect','Crystal seed']
    # skewness parameter --> a=0 is the normal distribution
    skewness = 12 # Skewness of the skewed Gaussian distribution
    fissure_region = (round(len(xv[0])/2)+2,100) # [0] middle point and [1] half width (nm)
    # Atomic specie -> Whay kind of atoms are affected by defects
    atomic_specie = 3 # Sulfur = 1 // Molibdenum = 3
    # Type of defects in the lattice
    defect_specie = 2 # Adatom = 2 // Crystal edge = 4 // Inner point of crystal = 5
    
    # 2 regions: etched and non-etched region
    # Boundary: 'vertical', 'horizontal', 'none', diagonal right, diagonal left
    mode = 3
    Boundary = ['vertical right', 'vertical left', 'horizontal', 'none','diagonal right', 'diagonal left']
    # Position: int - the row/column acting as a boundary and separe one region from the other
    # Position: coordinates when it is diagional boundary
    diagonal_right = [(int(i),int(i)) for i in np.arange(2 * device_size[0]/a - 1)]
    diagonal_left = [(int(2 * device_size[0]/a-1 - i),int(i)) for i in np.arange(2 * device_size[0]/a - 1)]
    Position = [round(len(xv[0])/2),round(len(xv[0])/2),round(len(xv[0])/2),round(len(xv[0])/2),diagonal_right,diagonal_left]
    split_regions = {'Boundary' : Boundary[mode], 'Position': Position[mode], 'ad_rate': parameters[1][n_sim]}
    

    prob_defects = parameters[0][n_sim]
    crystal_orientation = True
    pair_atom_defect=(3,4) # Introduce the crystal seed
    distribution_parameters = [distribution[5],skewness,fissure_region,pair_atom_defect,prob_defects,split_regions,rng]
    MoS2_lattice.defect_distributions(distribution_parameters)
    
    pair_atom_defect = (atomic_specie,defect_specie)
    distribution_parameters = [distribution[0],skewness,fissure_region,pair_atom_defect,prob_defects,split_regions,rng]

    
    MoS2_crystal = Cluster(MoS2_lattice.Grid_states) # Crystal seed
    
    MoS2_lattice.plot_lattice(crystal_orientation,'',0,0,True) # Initial state of the grid
    
    return MoS2_lattice,MoS2_crystal,distribution_parameters, paths,rng

def save_simulation(files_copy,dst,n_sim):
    

    parent_dir = 'Sim_'+str(n_sim)+'\\'
    os.makedirs(dst+parent_dir) 
    dst = dst+parent_dir
    program_directory = 'Program\\'
    data_directoy = 'Crystal growth\\'

    os.makedirs(dst + program_directory)
    os.makedirs(dst + data_directoy)
    
    paths = {'data': dst + data_directoy, 'program': dst + program_directory}

    for files in files_copy:
        shutil.copyfile(files, paths['program']+files)
        
    return paths

def save_variables(paths,variables):
    
    import shelve

    filename = 'variables'
    my_shelf = shelve.open(paths+filename,'n') # 'n' for new
    
    for key in variables:
        my_shelf[key] = variables[key]

    my_shelf.close()

