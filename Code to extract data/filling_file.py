# -*- coding: utf-8 -*-
"""
Created on Mon May 22 11:00:06 2023

@author: ALDANADS
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import pandas as pd
from hex_lattice import Hexagonal_lattice


name_generated_file = 'Figure.csv'

# Path to the variables files, for every folder and Sim_i
path = r'C:\Users\aldanads\OneDrive - TCDUD.onmicrosoft.com\2D device simulator project\Publications\Layer growth\Simulations\Sulfur limited process\20% of simulation domain\2 domains\Migration\Horizontal\\'
sampling_energy = True
if sampling_energy:
    folder_1 = os.listdir(path)
    
else:
    folder_1 = ['Vertical right_mig']
    
path_2 = r'\Program\\'

dfs = []
j = 0


for folder in folder_1:
    
    # Initialize variables
    Ad_rate = []
    Shaded_region = []
    Growth_t = []
    
    folders_sim = os.listdir(path+folder)
    
    for folder_sim in folders_sim:
    
        system = ['Windows','Linux']
        choose_system = system[1]
        
        if choose_system == 'Windows':
            import shelve
            
            filename = path + folder + r'\\' + folder_sim + path_2 + 'variables'

            my_shelf = shelve.open(filename)
            for key in my_shelf:
                globals()[key]=my_shelf[key]
            my_shelf.close()
            
        elif choose_system == 'Linux':
            
            import pickle
            filename = path + folder + r'\\' + folder_sim + path_2 + 'variables.pkl'
            print(filename)
            # Open the file in binary mode
            with open(filename, 'rb') as file:
              
                # Call load method to deserialze
                myvar = pickle.load(file)
                
            MoS2_crystal = myvar['MoS2_crystal']
            MoS2_lattice = myvar['MoS2_lattice']
        
        
        
        plt.rcParams["figure.dpi"] = 300
        
        n_Mo = MoS2_crystal.Mo_sites
        
        """
        Monolayers per second
        """
        # =============================================================================
        #print('Monolayer per second: ', sum(np.array(MoS2_lattice.v_adatom_flux)[0])/MoS2_lattice.time[-1]/n_Mo)
        #print('Adsorption-desorption equilibrium: ',sum(MoS2_lattice.v_adatom_flux[1])-MoS2_crystal.cluster_size)
        # 
        # =============================================================================
        Ad_rate.append(sum(np.array(MoS2_lattice.v_adatom_flux)[0])/MoS2_lattice.time[-1]/n_Mo)
        
        """
        Flake shape
        """
        flake_shape = False
        if flake_shape:
            aux = np.asarray(MoS2_crystal.cluster_ij)
            cluster_ij = ([aux[i][0] for i in np.arange(len(aux))],[aux[i][1] for i in np.arange(len(aux))])
            xv = MoS2_lattice.xv
            yv = MoS2_lattice.yv
            coord_xy_defects = MoS2_lattice.coord_xy_defects
            x_axis = MoS2_lattice.x_axis
            y_axis = MoS2_lattice.y_axis
            
            xv_1= xv[cluster_ij[0],cluster_ij[1]]
            yv_1= yv[cluster_ij[0],cluster_ij[1]]
            
            xv_2=xv[coord_xy_defects[0],coord_xy_defects[1]]
            yv_2= yv[coord_xy_defects[0],coord_xy_defects[1]]
                
            plt.scatter(xv[cluster_ij[0],cluster_ij[1]],yv[cluster_ij[0],cluster_ij[1]], color = 'black',s=5)
            plt.scatter(xv[coord_xy_defects[0],coord_xy_defects[1]],yv[coord_xy_defects[0],coord_xy_defects[1]], color = 'blue',s=5)
            plt.xlabel ("X axis (nm)")
            plt.ylabel ("Y axis (nm)")
            plt.xlim([-1,x_axis+1])
            plt.ylim([-1,y_axis+1])
            plt.show()
        """
        Proportion between etched and non-etched
        """
        proportion = True
        if proportion:
            cluster_ij = MoS2_crystal.cluster_ij
            xv = MoS2_lattice.xv
            
            position = round(len(xv[0])/2)
            
            # item = 0 is horizontal, item = 1 is vertical
            etched = [item for item in cluster_ij if item[0] >= position]
            
            area_etched = len(etched)
            
            # non_etched = [item for item in cluster_ij if item[0] < position]
            area_non_etched = len(cluster_ij) - len(etched)
            
            #print('Proportion of etched region: ',area_etched/len(cluster_ij))
            #print('Proportion of non-etched region: ',area_non_etched/len(cluster_ij))
            
            # Right
            Shaded_region.append(area_etched/len(cluster_ij))
            
            # Left
            #Shaded_region.append(1 - area_etched/len(cluster_ij))

        
        """
        Growth rate
        """
        growth_rate = True
        if growth_rate:
            Grid_states = MoS2_lattice.Grid_states
            x_axis = MoS2_lattice.x_axis
            y_axis = MoS2_lattice.y_axis
            area_Mo = x_axis*y_axis/n_Mo
            
            start = 1000
            
            fit_type = ['Exp','Linear']
            fitting = fit_type[0]
            
            if fitting == 'Exp':
            
                v_cluster_size = np.log(np.array(MoS2_crystal.v_cluster_size)*area_Mo)
                time = np.log(np.array(MoS2_lattice.time))
                
                model = LinearRegression().fit(time[start:].reshape(len(time)-start,1), v_cluster_size[start:])
                b0 = model.intercept_
                expb0 = np.exp(b0)
                b1 = model.coef_
                y_fit = expb0 * np.exp(b1 * time)
                
                plt.scatter(np.exp(time),np.exp(v_cluster_size), s = 10, label = 'Simulated data')
                plt.plot(np.exp(time),y_fit, color='black', linestyle='--', linewidth=2, label = 'Power law fitting')
                plt.legend()
                plt.xlabel('Time (s)')
                plt.ylabel('Flake area (nm$^2$)')
                plt.show()
                
                
                time_2 = np.log(np.arange(0,5*60,0.1))
                nm2_um2 = 1E-6
                y_fit_2 = expb0 * np.exp(b1 * time_2)
                plt.plot(np.exp(time_2),y_fit_2 * nm2_um2)
                plt.xlabel('Time (s)')
                plt.ylabel('Flake area (µm$^2$)')
                plt.show()
                
            elif fitting == 'Linear':
                
                v_cluster_size = np.array(MoS2_crystal.v_cluster_size)*area_Mo
                time = np.array(MoS2_lattice.time)
                
                model = LinearRegression().fit(time[start:].reshape(len(time)-start,1), v_cluster_size[start:])
                b0 = model.intercept_
                b1 = model.coef_
                y_fit =  b0 + b1*time
                
               
                plt.scatter(time,v_cluster_size, s = 10, label = 'Simulated data')
                plt.plot(time,y_fit, color='black', linestyle='--', linewidth=2, label = 'Power law fitting')
                plt.legend()
                plt.xlabel('Time (s)')
                plt.ylabel('Flake area (nm$^2$)')
                plt.show()
                
                
                time_2 = np.arange(0,5*60,0.1)
                nm2_um2 = 1E-6
                y_fit_2 = b0 + b1 * time_2
                
                plt.plot(time_2,y_fit_2 * nm2_um2)
                plt.xlabel('Time (s)')
                plt.ylabel('Flake area (µm$^2$)')
                plt.show()
                
                
            time = np.exp(time)
            #print(max(time))
            
            crystal_ij = np.where(Grid_states >= 4)
            xv = MoS2_lattice.xv
            yv = MoS2_lattice.yv
            coord_x = xv[crystal_ij[0],crystal_ij[1]]
            coord_y = yv[crystal_ij[0],crystal_ij[1]]
            
            Growth_t.append(max(time))

        
    #df = pd.DataFrame(Ad_rate,Shaded_region,Growth_t, columns=['Ad_rate_' + str(i), 'Shaded_region_' + str(i), 'Growth time_' + str(i)])
    
    df = pd.DataFrame({'Ad_rate_' + str(j+1): Ad_rate, 'Shaded_region_' + str(j+1): Shaded_region, 'Growth_time_' + str(j+1): Growth_t})   
    #df = pd.DataFrame({'Ad_rate_' + str(j+1): Ad_rate, 'Growth_time_' + str(j+1): Growth_t})
    j += 1
    dfs.append(df) # Data frame for each iteration i over each group of Sim
    
df_combined = pd.concat(dfs, axis=1)
    
df_combined.to_csv(name_generated_file, index=False)