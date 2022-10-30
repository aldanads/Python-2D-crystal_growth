# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 08:18:34 2022

@author: ALDANADS
"""
# ---------------------------------------------------------------------------#
########################### Hexagonal grid ###################################
# Reference: https://laurentperrinet.github.io/sciblog/posts/2020-04-16-creating-an-hexagonal-grid.html
# ---------------------------------------------------------------------------#

import numpy as np
import matplotlib.pyplot as plt

# We have a generic hexagonal grid --> We use this function to scale to the 
# MoS2 lattice size with lattice parameters
def scale_hex_grid(a,b,xv,yv):
    # We scale our grid to the size of the MoS2 lattice

    x_scale = np.sqrt(b**2-(a/2)**2)/1.5
    y_scale = a/2 
    xv=xv*x_scale
    yv=yv*y_scale
    
    # Check the angle between lattice constant --> It should be 60º for MoS2 
    # lattice
    #gamma=np.arctan(1.5*x_scale/y_scale)*180/np.pi
    #print (gamma)
    
    return xv,yv

# This function plot the Mo and the S in the hexagonal grid
def plot_Mo_S2(xv,yv,S_color,Mo_color):
    
    # Sulfur atomic radius: 100 pm
    # Moldibdenum atomic radius: 139 pm
    # Change s to make them in scale
    plt.scatter(xv[1::2,0::3],yv[1::2,0::3],s=5,color = S_color,) # Sulfur
    plt.scatter(xv[0::2,1::3],yv[0::2,1::3],s=5,color = S_color,) # Sulfur
    plt.scatter(xv[0::2,0::3],yv[0::2,0::3],s=1.39,color = Mo_color,) # Molibdenum
    plt.scatter(xv[1::2,2::3],yv[1::2,2::3],s=1.39,color = Mo_color,) # Molibdenum
    
    return
    
# ------------------------------- Rectangular grid ---------------------------#
"""
N=121
N_X = int(np.sqrt(N))
#Floor Division (also called Integer Division): rounded to the next smallest whole number
N_Y = N // N_X
xv, yv = np.meshgrid(np.arange(N_X), np.arange(N_Y), sparse = False, indexing='xy')
#fig,ax = plt.subplots(figsize=(fig_width,fig_width)
plt.scatter(xv,yv)
plt.show()
"""
# ----------------------------------------------------------------------------#

# ------------------------------- Hexagonal grid -----------------------------#

# Create hegaxonal grid
def create_hex_grid(a,b,device_size,plot_grid):
    """
    Mortazavi, Majid, Chao Wang, Junkai Deng, Vivek B. Shenoy, and Nikhil V. Medhekar. 
    "Ab initio characterization of layered MoS2 as anode for sodium-ion batteries." 
    Journal of Power Sources 268 (2014): 279-286.
    """
    
    S_color='orange'
    Mo_color='purple'
    #  axes are proportional to a/2 --> Multiply by 2 to have nm
    grid_size = (device_size[0] * 2, device_size[1] * 2) 
    N_X=int(grid_size[0]/a)
    N_Y=int(grid_size[1]/a)

    # Lattice size
    #N=121
    ratio=np.sqrt(3)/2 
    #N_X = int(np.sqrt(N)/ratio)
    #N_Y = N // N_X
    #N_X=grid_size(0)/a
    xv, yv = np.meshgrid(np.arange(N_X), np.arange(N_Y), sparse=False, indexing='xy')
    xv = xv * ratio
    xv[::2] += ratio/2

    xv, yv = scale_hex_grid(a,b,xv,yv)
    
    # Plot the grid we created
    if plot_grid == True:
        plot_Mo_S2(xv,yv,S_color,Mo_color)

        # Generic hexagonal grid
        #plt.scatter(xv,yv,s=1)
        plt.xlabel ("X axis (nm)")
        plt.ylabel ("Y axis (nm)")

        plt.show()
        
    return xv,yv






