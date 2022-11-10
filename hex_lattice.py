# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 18:22:27 2022

@author: ALDANADS
"""

# Create hegaxonal grid --> Object

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import skewnorm



class Hexagonal_lattice():
    
    def __init__(self,a,b,device_size,atom_colors,Act_E,T):
        self.a = a
        self.b = b
        self.x_axis = device_size[0]
        self.y_axis = device_size[1]
        self.atom_colors = atom_colors
        self.Act_E = Act_E
        self.T = T
        
        self.create_hex_grid()
        self.pristine_crystal()

    
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
            
        
    
    def plot_lattice(self,crystal_orientation = False):
        
        # Sulfur atomic radius: 100 pm
        # Moldibdenum atomic radius: 139 pm
        # Change s to make them in scale
        plt.scatter(self.xv[1::2,0::3],self.yv[1::2,0::3],color = self.atom_colors[0],s=1) # Sulfur
        plt.scatter(self.xv[0::2,1::3],self.yv[0::2,1::3],color = self.atom_colors[0],s=1) # Sulfur
        plt.scatter(self.xv[0::2,0::3],self.yv[0::2,0::3],color = self.atom_colors[1],s=1.39) # Molibdenum
        plt.scatter(self.xv[1::2,2::3],self.yv[1::2,2::3],color = self.atom_colors[1],s=1.39) # Molibdenum
        
        if (type(self.Grid_states) == np.ndarray):
            coord_xy_Vs = np.where(self.Grid_states == 2)
            plt.scatter(self.xv[coord_xy_Vs[0],coord_xy_Vs[1]],self.yv[coord_xy_Vs[0],coord_xy_Vs[1]], color = self.atom_colors[2],s=5)
            coord_xy_Vs = np.where(self.Grid_states == 4)
            plt.scatter(self.xv[coord_xy_Vs[0],coord_xy_Vs[1]],self.yv[coord_xy_Vs[0],coord_xy_Vs[1]], color = self.atom_colors[3],s=5)

        if (crystal_orientation == True): # Crystal orientation
            arrow1 = plt.arrow(self.xv[2,0],self.yv[2,0],self.x_axis/4,0,width =0.05)
            arrow2 = plt.arrow(self.xv[2,0],self.yv[2,0],0,self.y_axis/4,width =0.05, color = 'green')
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
        
        
        """
        -------------------- Introducing defects in the crystal ---------------
        """
        
    def pristine_crystal(self):
        
        Grid_states=np.zeros((len(self.xv),len(self.xv[0])), dtype=int)
        
        Grid_states[1::2,0::3]=1 # Sulfur
        Grid_states[0::2,1::3]=1 # Sulfur
        Grid_states[0::2,0::3]=3 # Molibdenum
        Grid_states[1::2,2::3]=3 # Molibdenum
        
        self.Grid_states = Grid_states
        
    def introduce_defects_j_row(self,j,prob_defects):
        
        counter=0

        # The defects we are introducing in column j
        prob_defects=np.random.rand(sum(self.Grid_states[:,j] == self.atomic_specie)) < prob_defects
            
        for i in np.arange(len(self.xv)):
            if self.Grid_states[i,j] == self.atomic_specie:
                
                if prob_defects[counter]:
                    self.Grid_states[i,j] = 2
                
                counter += 1

    # Uniform distribution
    def defects_row(self,prob_defects,fissure_region):
        
        
        # 1 position is sulfur and the other is Molybdenum --> there is len(xv)/2 sulfur in a column
        
        a = (self.xv[0,1]-self.xv[0,0])*2 # lattice constant a (nm)

        irradiated_row = fissure_region[0] # Middle point of the triangle base
        # Half of the triangle base
        width_fissure = fissure_region[1]*2/a # half width of fissure region in columns
        
        length_xv = len(self.xv[0])
        # If the fissure region is greater than the simulation domain, we cover the simulation domain
        if 2 * width_fissure > length_xv:
            start_row = 0
            finish_row = length_xv
        else: # The fissure region fit the simulation domain
            # The triangle starting point in the grid
            start_row = round(irradiated_row - width_fissure)
            # The last point of the triangle in the grid
            finish_row = round(irradiated_row + width_fissure)
            
        list_prob = np.zeros(length_xv)
        list_prob[start_row:finish_row] = prob_defects
        
        self.list_prob = list_prob
        
        for j in np.arange(start_row,finish_row):
            
            self.introduce_defects_j_row(j,prob_defects)

    # Skewed Gaussian distribution:
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.skewnorm.html
    def defects_skewed_gaussian(self,prob_defects,fissure_region,skewness):

        a = (self.xv[0,1]-self.xv[0,0])*2 # lattice constant a (nm)
        
        irradiated_row = fissure_region[0]
        width_fissure = fissure_region[1]*2/a # width of fissure region in columns
        length_xv = len(self.xv[0])

        
        mean = skewnorm.mean(skewness, moments='m') # mean of the distribution
        std = skewnorm.std(skewness) # standard deviation
        x = np.linspace(skewnorm.ppf(0.001, skewness), skewnorm.ppf(0.999, skewness), length_xv)
        norm_const=sum(skewnorm.pdf(x, skewness)) # Normalization constant
        prob = skewnorm.pdf(x, skewness)/norm_const # Probability density function normalized to 1
        max_prob = max(prob) # Maximum probability found in the distribution
        max_prob_defect = prob_defects/max_prob # Scale the distribution
        prob = max_prob_defect*prob # We set the peak at the max probability of generating a vacancy
        
        list_prob = []
        
        for j in np.arange(0,length_xv):
            
            # Location of the peak density in the simulation domain
            x=mean-(irradiated_row-j)*std/width_fissure;
            
            prob = max_prob_defect * skewnorm.pdf(x, skewness) / norm_const
            
            list_prob.append(prob)

            self.introduce_defects_j_row(j,prob_defects)
            
        self.list_prob = list_prob
        
    def defect_triangle(self,prob_defects,fissure_region):
                    
        a = (self.xv[0,1]-self.xv[0,0])*2 # lattice constant a (nm)
             
        irradiated_row = fissure_region[0] # Middle point of the triangle base
        # Half of the triangle base
        width_fissure = fissure_region[1]*2/a # half width of fissure region in columns
            
        # The triangle starting point in the grid
        start_triangle = round(irradiated_row - width_fissure)
        # The last point of the triangle in the grid
        finish_triangle = round(irradiated_row + width_fissure)
            
        # Slope --> right triangle hypotenuse 
        slope = prob_defects / (2 * width_fissure)
        b = prob_defects
            
        list_prob = np.zeros(len(self.xv[0]))
            
        for j in np.arange(start_triangle,finish_triangle):
            # Triangle hypotenuse 
            dose=slope*(start_triangle-j)+b;
            list_prob[j] = dose

            self.introduce_defects_j_row(j,prob_defects)
            
        self.list_prob = list_prob
        
    def defect_distributions(self,prob_defects,fissure_region,skewness,distribution,pair_atom_defect):
        
        self.atomic_specie = pair_atom_defect[0]
        self.defect_specie = pair_atom_defect[1]
        
        if distribution == 'uniform':
            self.defects_row(prob_defects,fissure_region)
        if distribution == 'triangle':
            self.defect_triangle(prob_defects,fissure_region)
        if distribution == 'skewed_gaussian':
            self.defects_skewed_gaussian(prob_defects,fissure_region,skewness)
        if distribution == 'test 1: single adatom':
            self.single_defect()
            
    def single_defect(self):

        j = 0
        length_xv = len(self.xv)
        length_yv = len(self.yv[0])

        if self.Grid_states[int(length_xv/2),int(length_yv/2)] == self.atomic_specie:
            self.Grid_states[int(length_xv/2),int(length_yv/2)] = self.defect_specie
        else:
            while self.Grid_states[int(length_xv/2),int(length_yv/2)+j] != self.atomic_specie:
                j += 1
            self.Grid_states[int(length_xv/2),int(length_yv/2)+j] = self.defect_specie 
            self.Grid_states[int(length_xv/2)+1,int(length_yv/2)+j+1] = self.defect_specie 

    
            
    def coord_defects(self):
        
        self.n_defects = np.count_nonzero(self.Grid_states == 2)
        self.coord_xy_defects = np.where(self.Grid_states == 2)
        
        
        
        
        
            

