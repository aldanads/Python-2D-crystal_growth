# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 20:42:09 2022

@author: ALDANADS
"""
import numpy as np
from scipy.stats import skewnorm

# Initialize the sulfurs and molibdenums in the lattice
def defects(xv,yv):
    
    Grid_states=np.zeros((len(xv),len(xv[0])), dtype=int)
    
    Grid_states[1::2,0::3]=1 # Sulfur
    Grid_states[0::2,1::3]=1 # Sulfur
    Grid_states[0::2,0::3]=3 # Molibdenum
    Grid_states[1::2,2::3]=3 # Molibdenum
    
    return Grid_states
    
# Uniform distribution of defects    
def defects_in_row(xv,yv,prob_gen_defect,fissure_region):
    
    Grid_states = defects(xv,yv) # Initialize a pristine 2D MoS2 layer
    # 1 position is sulfur and the other is Molybdenum --> there is len(xv)/2 sulfur in a column
    
    a = (xv[0,1]-xv[0,0])*2 # lattice constant a (nm)

    irradiated_row = fissure_region[0] # Middle point of the triangle base
    # Half of the triangle base
    width_fissure = fissure_region[1]*2/a # half width of fissure region in columns
    
    # If the fissure region is greater than the simulation domain, we cover the simulation domain
    if 2 * width_fissure > len(xv[0]):
        start_row = 0
        finish_row = len(xv[0])
    else: # The fissure region fit the simulation domain
        # The triangle starting point in the grid
        start_row = round(irradiated_row - width_fissure)
        # The last point of the triangle in the grid
        finish_row = round(irradiated_row + width_fissure)
        
    list_dose = np.zeros(len(xv[0]))
    list_dose[start_row:finish_row] = prob_gen_defect
    
    for j in np.arange(start_row,finish_row):
        
        Grid_states = introduce_defects_j_row(Grid_states,j,prob_gen_defect,xv)
            
    return Grid_states, list_dose

# Skewed Gaussian distribution:
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.skewnorm.html
def defects_skewed_gaussian(xv,yv,prob_gen_defect,fissure_region,skewness):
    
    Grid_states = defects(xv,yv) # Initialize a pristine 2D MoS2 layer

    a = (xv[0,1]-xv[0,0])*2 # lattice constant a (nm)
    
    irradiated_row = fissure_region[0]
    width_fissure = fissure_region[1]*2/a # width of fissure region in columns
    
    mean = skewnorm.mean(skewness, moments='m') # mean of the distribution
    std = skewnorm.std(skewness) # standard deviation
    x = np.linspace(skewnorm.ppf(0.001, skewness), skewnorm.ppf(0.999, skewness), len(xv[0]))
    norm_const=sum(skewnorm.pdf(x, skewness)) # Normalization constant
    dose = skewnorm.pdf(x, skewness)/norm_const # Probability density function normalized to 1
    max_prob = max(dose) # Maximum probability found in the distribution
    max_dose = prob_gen_defect/max_prob # Scale the distribution
    dose = max_dose*dose # We set the peak at the max probability of generating a vacancy
    
    list_dose = []
    
    for j in np.arange(0,len(xv[0])):
        
        # Location of the peak density in the simulation domain
        x=mean-(irradiated_row-j)*std/width_fissure;
        
        dose = max_dose * skewnorm.pdf(x, skewness) / norm_const
        
        list_dose.append(dose)

        Grid_states = introduce_defects_j_row(Grid_states,j,dose,xv)

    
    return Grid_states,list_dose

# Right triangle distribution of defects:
def defect_triangle(xv,yv,prob_gen_defect,fissure_region):
    
    Grid_states = defects(xv,yv) # Initialize a pristine 2D MoS2 layer
    a = (xv[0,1]-xv[0,0])*2 # lattice constant a (nm)
     
    irradiated_row = fissure_region[0] # Middle point of the triangle base
    # Half of the triangle base
    width_fissure = fissure_region[1]*2/a # half width of fissure region in columns
    
    # The triangle starting point in the grid
    start_triangle = round(irradiated_row - width_fissure)
    # The last point of the triangle in the grid
    finish_triangle = round(irradiated_row + width_fissure)
    
    # Slope --> right triangle hypotenuse 
    slope = prob_gen_defect / (2 * width_fissure)
    b = prob_gen_defect
    
    list_dose = np.zeros(len(xv[0]))
    
    for j in np.arange(start_triangle,finish_triangle):
        # Triangle hypotenuse 
        dose=slope*(start_triangle-j)+b;
        list_dose[j] = dose

        Grid_states = introduce_defects_j_row(Grid_states,j,dose,xv)

    
    return Grid_states, list_dose

# Selecting the kind of distribution
def defect_distributions(xv,yv,prob_gen_defect,fissure_region,skewness,distribution):
    
    if distribution == 'uniform':
        Grid_states,list_dose = defects_in_row(xv,yv,prob_gen_defect,fissure_region)
    if distribution == 'triangle':
        Grid_states,list_dose = defect_triangle(xv,yv,prob_gen_defect,fissure_region)
    if distribution == 'skewed_gaussian':
        Grid_states,list_dose = defects_skewed_gaussian(xv,yv,prob_gen_defect,fissure_region,skewness)

    return Grid_states, list_dose

# Introduce defects in row j dependending on the probability dose
def introduce_defects_j_row(Grid_states,j,dose,xv):
    
    counter=0

    # The defects we are introducing in column j
    prob_defects=np.random.rand(sum(Grid_states[:,j] == 1)) < dose
        
    for i in np.arange(len(xv)):
        if Grid_states[i,j] == 1:
            
            if prob_defects[counter]:
                Grid_states[i,j] = 2
            
            counter += 1
    
    return Grid_states


