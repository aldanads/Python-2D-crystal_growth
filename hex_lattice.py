# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 18:22:27 2022

@author: ALDANADS
"""

# Create hegaxonal grid --> Object

import numpy as np
import matplotlib.pyplot as plt


class Hexagonal_lattice():
    
    def __init__(self,a,b,device_size,atom_colors):
        self.a = a
        self.b = b
        self.x_axis = device_size[0]
        self.y_axis = device_size[1]
        self.atom1_colors = atom_colors[0]
        self.atom2_colors = atom_colors[1]
        self.atom3_colors = atom_colors[2]

    
    def create_hex_grid(self):

        """
        Mortazavi, Majid, Chao Wang, Junkai Deng, Vivek B. Shenoy, and Nikhil V. Medhekar. 
        "Ab initio characterization of layered MoS2 as anode for sodium-ion batteries." 
        Journal of Power Sources 268 (2014): 279-286.
        """
        
        #  axes are proportional to a/2 --> Multiply by 2 to have nm
        grid_size = (self.x_axis * 2, self.y_axis * 2) 
        N_X=int(grid_size[0]/self.a)
        N_Y=int(grid_size[1]/self.a)

        # Lattice size
        #N=121
        ratio=np.sqrt(3)/2 
        #N_X = int(np.sqrt(N)/ratio)
        #N_Y = N // N_X
        #N_X=grid_size(0)/a
        xv, yv = np.meshgrid(np.arange(N_X), np.arange(N_Y), sparse=False, indexing='xy')
        xv = xv * ratio
        xv[::2] += ratio/2

        # We scale our grid to the size of the lattice defined by a and b
        x_scale = np.sqrt(self.b**2-(self.a/2)**2)/1.5
        y_scale = self.a/2 
        xv=xv*x_scale
        yv=yv*y_scale
        
        self.xv = xv
        self.yv = yv
            
        
    
    def plot_lattice(self,Grid_states = 0,crystal_orientation = False):
        
        # Sulfur atomic radius: 100 pm
        # Moldibdenum atomic radius: 139 pm
        # Change s to make them in scale
        plt.scatter(self.xv[1::2,0::3],self.yv[1::2,0::3],color = self.atom1_colors,s=1) # Sulfur
        plt.scatter(self.xv[0::2,1::3],self.yv[0::2,1::3],color = self.atom1_colors,s=1) # Sulfur
        plt.scatter(self.xv[0::2,0::3],self.yv[0::2,0::3],color = self.atom2_colors,s=1.39) # Molibdenum
        plt.scatter(self.xv[1::2,2::3],self.yv[1::2,2::3],color = self.atom2_colors,s=1.39) # Molibdenum
        
        if (type(Grid_states) == np.ndarray):
            coord_xy_Vs = np.where(Grid_states == 2)
            plt.scatter(self.xv[coord_xy_Vs[0],coord_xy_Vs[1]],self.yv[coord_xy_Vs[0],coord_xy_Vs[1]], color = self.atom3_colors,s=1)
        if (crystal_orientation == True):
            arrow1 = plt.arrow(self.xv[2,0],self.yv[2,0],1,0,width =0.05)
            arrow2 = plt.arrow(self.xv[2,0],self.yv[2,0],0,1,width =0.05, color = 'green')
            plt.legend([arrow1,arrow2], ['Armchair','Zigzag'])
        plt.xlabel ("X axis (nm)")
        plt.ylabel ("Y axis (nm)")

        plt.show()
    
    
    # Generic hexagonal grid
    def plot_generic_hex_grid(self):
        plt.scatter(self.xv,self.yv,s=1)
        plt.xlabel ("X axis (nm)")
        plt.ylabel ("Y axis (nm)")

        plt.show()


