# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 20:15:19 2022

@author: ALDANADS
"""

from hex_grid import*
from defects import*
import numpy as np
import matplotlib.pyplot as plt


""" --------------------------------------------------------------------------
 -----------------------Creating of the hexagonal grid -----------------------
 -----------------------------------------------------------------------------
""" 
a=0.639 # nm
b=0.639 # nm
plot_grid = True # Boolean to plot the hex grid we created
device_size = (50 , 50) # Size of the grid in nm. grid_size[0] is x and grid_size[1] is y.


# Creation of hexagonal grid
# xv and yv is the coordinates of the different points of the hexagonal grid
xv,yv = create_hex_grid(a,b,device_size,plot_grid) 

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
