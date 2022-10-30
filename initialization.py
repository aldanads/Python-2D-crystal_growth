# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 20:15:19 2022

@author: ALDANADS
"""

from hex_grid import*
from defects import*
from hex_lattice import*
import numpy as np
import matplotlib.pyplot as plt


""" --------------------------------------------------------------------------
 -----------------------Creating of the hexagonal grid -----------------------
 -----------------------------------------------------------------------------
""" 
a=0.639 # nm
b=0.639 # nm
device_size = (15 , 15) # Size of the grid in nm. grid_size[0] is x and grid_size[1] is y.
atom_colors=['orange','purple'] # MoS2 -> First is Sulfur and second is Mo


MoS2_lattice = Hexagonal_lattice(a,b,device_size,atom_colors)
MoS2_lattice.create_hex_grid()
xv=MoS2_lattice.xv
yv=MoS2_lattice.yv

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

"""
distribution:
    'uniform' --> Uniform rows
    'skewed_gaussian' --> Skewed Gaussian distribution of defects
    'triangle' --> Right triangle distribution of defects
"""


distribution = 'uniform'
prob_gen_defect = 0.5 # prob of generating defects --> Peak density
fissure_region = (round(len(xv[0])/2)+1,4) # [0] middle point and [1] half width (nm)

# skewness parameter --> a=0 is the normal distribution
skewness = 12 # Skewness of the skewed Gaussian distribution

Grid_states, list_dose = defect_distributions(xv,yv,prob_gen_defect,fissure_region,skewness,distribution)
