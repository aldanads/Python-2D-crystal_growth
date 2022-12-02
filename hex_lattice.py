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
    
    def __init__(self,a,b,device_size,atom_colors,Backup_energy,T):
        self.a = a
        self.b = b
        self.x_axis = device_size[0]
        self.y_axis = device_size[1]
        self.atom_colors = atom_colors
        self.Backup_energy = Backup_energy
        self.T = T
        self.time = [0]
        
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
            
        
    
    def plot_lattice(self,crystal_orientation = False, path = '',t=0,i=0, grid=False,cluster_ij = False,split_regions = False):
        
        # Sulfur atomic radius: 100 pm
        # Moldibdenum atomic radius: 139 pm
        # Change s to make them in scale
        
            if (grid == True) or self.x_axis <= 50:
                plt.scatter(self.xv[1::2,0::3],self.yv[1::2,0::3],color = self.atom_colors[0],s=1) # Sulfur
                plt.scatter(self.xv[0::2,1::3],self.yv[0::2,1::3],color = self.atom_colors[0],s=1) # Sulfur
                plt.scatter(self.xv[0::2,0::3],self.yv[0::2,0::3],color = self.atom_colors[1],s=1.39) # Molibdenum
                plt.scatter(self.xv[1::2,2::3],self.yv[1::2,2::3],color = self.atom_colors[1],s=1.39) # Molibdenum
                
            if split_regions != False:
                
                split_line = split_regions['Position'] # The boundary position
                if split_regions['Boundary'] == 'vertical':
                    y = [self.yv[len(self.yv)-1,0]] * (len(self.yv)-split_line) # max value of y
                    # x go from the split line until the end
                    plt.fill_between(self.xv[0,split_line:],y,alpha = 0.2, color = 'blue')
                elif split_regions['Boundary'] == 'horizontal':
                    x = [self.xv[0,len(self.xv)-1]] * (len(self.xv)-split_line) # max value of x
                    # y go from the split line until the end
                    plt.fill_betweenx(self.yv[split_line:,0],x,alpha = 0.2, color = 'blue')
                
            
            if (type(self.Grid_states) == np.ndarray):
                
                if hasattr(self,'coord_xy_defects') == False:
                    coord_xy_defects = np.where(self.Grid_states == 2)
                    plt.scatter(self.xv[coord_xy_defects[0],coord_xy_defects[1]],self.yv[coord_xy_defects[0],coord_xy_defects[1]], color = self.atom_colors[2],s=5)
                else:
                    plt.scatter(self.xv[self.coord_xy_defects[0],self.coord_xy_defects[1]],self.yv[self.coord_xy_defects[0],self.coord_xy_defects[1]], color = self.atom_colors[2],s=5)
                
                
                if cluster_ij == False: cluster_ij = np.where(self.Grid_states >= 4)
                else: 
                    # Transform the list of tuples to a tuple with two lists
                    aux = np.asarray(cluster_ij)
                    cluster_ij = ([aux[i][0] for i in np.arange(len(aux))],[aux[i][1] for i in np.arange(len(aux))])
                
                #print(cluster_ij)

                plt.scatter(self.xv[cluster_ij[0],cluster_ij[1]],self.yv[cluster_ij[0],cluster_ij[1]], color = self.atom_colors[3],s=5)

            if (crystal_orientation == True): # Crystal orientation
                arrow1 = plt.arrow(self.xv[2,0],self.yv[2,0],self.x_axis/4,0,width =0.05)
                arrow2 = plt.arrow(self.xv[2,0],self.yv[2,0],0,self.y_axis/4,width =0.05, color = 'green')
                plt.legend([arrow1,arrow2], ['Armchair','Zigzag'])
            plt.xlabel ("X axis (nm)")
            plt.ylabel ("Y axis (nm)")
            
            if path == '':
                plt.show()
            else:
                plt.savefig(path+str(i)+'_t(s) = '+str(round(t,5))+' .png', dpi = 100)
    
    
    # Generic hexagonal grid
    def plot_generic_hex_grid(self):
        plt.scatter(self.xv,self.yv,s=1)
        plt.xlabel ("X axis (nm)")
        plt.ylabel ("Y axis (nm)")

        plt.show()
        
        
        """
        -------------------- Introducing defects in the crystal ---------------
        """
        
    def defect_distributions(self,distribution_parameters):
            
        distribution = distribution_parameters[0]
        pair_atom_defect = distribution_parameters[3]
        rng = distribution_parameters[6]

        
        self.atomic_specie = pair_atom_defect[0]
        self.defect_specie = pair_atom_defect[1]
            
        if distribution == 'uniform':
            fissure_region = distribution_parameters[2] 
            prob_defects = distribution_parameters[4]
            split_regions = distribution_parameters[5]
            self.defects_row(prob_defects,fissure_region,rng,split_regions)
            
        if distribution == 'triangle':
            fissure_region = distribution_parameters[2] 
            prob_defects = distribution_parameters[4]
            self.defect_triangle(prob_defects,fissure_region,rng)
            
        if distribution == 'skewed_gaussian':
            fissure_region = distribution_parameters[2] 
            skewness = distribution_parameters[1] 
            prob_defects = distribution_parameters[4]
            self.defects_skewed_gaussian(prob_defects,fissure_region,skewness,rng)
            
        if distribution == 'test 1: single adatom':
            self.single_defect(rng)
            
        if distribution == 'test 2: column defect':
            self.column_defect(rng)
            
        if distribution == 'Crystal seed':
            self.crystal_seed()
            
        
    def pristine_crystal(self):
        
        Grid_states=np.zeros((len(self.xv),len(self.xv[0])), dtype=int)
        
        Grid_states[1::2,0::3]=1 # Sulfur
        Grid_states[0::2,1::3]=1 # Sulfur
        Grid_states[0::2,0::3]=3 # Molibdenum
        Grid_states[1::2,2::3]=3 # Molibdenum
        
        self.Grid_states = Grid_states
        

    # Uniform distribution
    def defects_row(self,prob_defects,fissure_region,rng,split_regions):
        
        
        # 1 position is sulfur and the other is Molybdenum --> there is len(xv)/2 sulfur in a column
        
        a = (self.xv[0,1]-self.xv[0,0])*2 # lattice constant a (nm)

        irradiated_row = fissure_region[0] # Middle point of the triangle base
        # Half of the triangle base
        width_fissure = fissure_region[1]*2/a # half width of fissure region in columns
        
        length_xv = len(self.xv[0])
        self.adatom_flux = 0
        list_prob = np.zeros(length_xv)

        if split_regions['Boundary'] == 'none': # Only one region
            # If the fissure region is greater than the simulation domain, we cover the simulation domain
            if 2 * width_fissure > length_xv:
                start_row = 0
                finish_row = length_xv
            else: # The fissure region fit the simulation domain
                # The triangle starting point in the grid
                start_row = round(irradiated_row - width_fissure)
                # The last point of the triangle in the grid
                finish_row = round(irradiated_row + width_fissure)
                
                list_prob[start_row:finish_row] = prob_defects
            
            for j in np.arange(start_row,finish_row):
                
                self.introduce_defects_j_row(j,prob_defects,rng)
                
        elif split_regions['Boundary'] == 'vertical': # Two regions splitted by a vertical boundary
            start_row = 0
            finish_row = length_xv    
            list_prob[0:split_regions['Position']] = prob_defects
            list_prob[split_regions['Position']:] = split_regions['ad_rate']
            
            for j in np.arange(start_row,finish_row):
                
                self.introduce_defects_j_row(j,list_prob[j],rng)
            
        elif split_regions['Boundary'] == 'horizontal': # Two regions splitted by a horizontal boundary
            start_row = 0
            finish_row = length_xv  
            
            for j in np.arange(start_row,finish_row):
                
                self.introduce_defects_j_row(j,prob_defects,rng,split_regions)

        
        self.list_prob = list_prob
        

            
        # adatom flux --> adatoms adsorted by the substrate (adatoms / nm^2)
        self.adatom_flux = self.adatom_flux / (self.x_axis * self.y_axis)
        
    # Triangle distribution of defects
    def defect_triangle(self,prob_defects,fissure_region,rng):
                    
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

            self.introduce_defects_j_row(j,prob_defects,rng)
            
        self.list_prob = list_prob

    # Skewed Gaussian distribution:
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.skewnorm.html
    def defects_skewed_gaussian(self,prob_defects,fissure_region,skewness,rng):

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

            self.introduce_defects_j_row(j,prob_defects,rng)
            
        self.list_prob = list_prob
        
            
    def single_defect(self,rng):

        j = 0
        length_xv = len(self.xv)
        length_yv = len(self.yv[0])
        x = int(length_xv/2)
        y = int(length_yv/2)

        if self.Grid_states[x,y] == self.atomic_specie:
            self.Grid_states[x,y] = self.defect_specie

        else:
            while self.Grid_states[x,y+j] != self.atomic_specie:
                j += 1
                
            self.Grid_states[x,y+j] = self.defect_specie 

            
    def column_defect(self,rng):
        
        j = 0
        length_xv = len(self.xv)
        length_yv = len(self.yv[0])
        x = int(length_xv/2)
        y = int(length_yv/2)

        if self.Grid_states[x,y] == self.atomic_specie:
            self.Grid_states[x,y] = self.defect_specie
            self.introduce_defects_j_row(y+j,1,rng)
            self.introduce_defects_j_row(y+j-1,1,rng)
            self.Grid_states[x,y+j-3] = 4

        else:
            while self.Grid_states[x,y+j] != self.atomic_specie:
                j += 1
                
            #self.Grid_states[x,y+j] = self.defect_specie 
            self.introduce_defects_j_row(y+j,1,rng)
            self.introduce_defects_j_row(y+j+1,1,rng)
            self.Grid_states[x-1,y+j-2] = 4

    def crystal_seed(self):
        j = 0
        length_xv = len(self.xv)
        length_yv = len(self.yv[0])
        x = int(length_xv/2)
        y = int(length_yv/2)
        
        if self.Grid_states[x,y] == self.atomic_specie:
            if self.Grid_states[x-1,y] == 0:
                self.Grid_states[x,y] = self.defect_specie 
                self.Grid_states[x+1,y+1] = self.defect_specie
                
                self.Grid_states[x+2,y-3] = self.defect_specie
                self.Grid_states[x-2,y-3] = self.defect_specie
                self.Grid_states[x,y-3] = self.defect_specie 

                self.Grid_states[x-1,y-2] = self.defect_specie
                self.Grid_states[x+1,y-2] = self.defect_specie
            else:
                self.Grid_states[x,y] = self.defect_specie 
                self.Grid_states[x+1,y+2] = self.defect_specie
                
                self.Grid_states[x+2,y-3] = self.defect_specie
                self.Grid_states[x-2,y-3] = self.defect_specie
                self.Grid_states[x,y-3] = self.defect_specie 

                self.Grid_states[x-1,y-1] = self.defect_specie
                self.Grid_states[x+1,y-1] = self.defect_specie
                

        else:
            while self.Grid_states[x,y+j] != self.atomic_specie:
                j += 1
                

            if self.Grid_states[x-1,y+j] == 0:
  
                self.Grid_states[x,y+j] = self.defect_specie 
                self.Grid_states[x+1,y+j+1] = self.defect_specie
                
                self.Grid_states[x+2,y+j-3] = self.defect_specie
                self.Grid_states[x-2,y+j-3] = self.defect_specie
                self.Grid_states[x,y+j-3] = self.defect_specie 

                self.Grid_states[x-1,y+j-2] = self.defect_specie
                self.Grid_states[x+1,y+j-2] = self.defect_specie
                
            else:
                
                self.Grid_states[x,y+j] = self.defect_specie 
                self.Grid_states[x+1,y+j+2] = self.defect_specie
                self.Grid_states[x+1,y+j+1] = self.defect_specie

                
                self.Grid_states[x+2,y+j-3] = self.defect_specie
                self.Grid_states[x-2,y+j-3] = self.defect_specie
                self.Grid_states[x,y-3] = self.defect_specie 

                self.Grid_states[x-1,y+j-1] = self.defect_specie
                self.Grid_states[x+1,y+j-1] = self.defect_specie
            
        
            
    def coord_defects(self):
        
        n_defects = np.count_nonzero(self.Grid_states == 2)
        coord_xy_defects = np.where(self.Grid_states == 2)
        self.coord_xy_defects = coord_xy_defects
        self.n_defects = n_defects
        # Transform the coordinate in a list of tuples --> Easier to compare with other tuples
        coord_xy_defects = [(coord_xy_defects[0][i],coord_xy_defects[1][i]) for i in np.arange(len(coord_xy_defects[0]))]
        
        return coord_xy_defects
        
    def introduce_defects_j_row(self,j,prob_defects,rng, split_regions = False):
        
        counter=0
        x_length = len(self.xv)
        #init_x = x_length/2*
        # The defects we are introducing in column j
        if split_regions == False:
            defects_j_column = rng.random(sum(self.Grid_states[:,j] == self.atomic_specie)) < prob_defects
        elif split_regions['Boundary'] == 'horizontal':  
            # Two lists: Before the split line and after the split line
            defects_j_column = list(rng.random(sum(self.Grid_states[:split_regions['Position'],j] == self.atomic_specie)) < prob_defects) 
            defects_j_column_aux = list(rng.random(sum(self.Grid_states[split_regions['Position']:,j] == self.atomic_specie)) < split_regions['ad_rate'])
            defects_j_column += defects_j_column_aux
            
        self.adatom_flux += sum(defects_j_column)
        
        for i in np.arange(x_length):
            if self.Grid_states[i,j] == self.atomic_specie:
                
                if defects_j_column[counter]:
                    self.Grid_states[i,j] = self.defect_specie
                
                counter += 1
   
    def add_time(self,t):
        
        self.time.append(self.time[len(self.time)-1] + t)
        
            
            
    


        
        
        
        
        
            

