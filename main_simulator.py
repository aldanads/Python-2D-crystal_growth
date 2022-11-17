# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:01:13 2022

@author: ALDANADS
"""
from initialization import*
from KMC import*
import shutil

E_propagation = [1.4,1.4,1.4,1.4,1.4,1.3,1.4,1.5,1.6,1.7]
E_nucleation = [1.7,1.6,1.5,1.4,1.3,1.7,1.7,1.7,1.7,1.7]

E_nuc_prop = [E_nucleation,E_propagation]

for n_sim in np.arange(len(E_nuc_prop[0])):

    MoS2_lattice, MoS2_crystal,distribution_parameters,dst_data = initialization(E_nuc_prop,n_sim)
    prob = np.zeros(len(MoS2_lattice.Act_E)+1)
    
    for i in np.arange(12000):
        MoS2_lattice, MoS2_crystal,Mo_adatom,prob = KMC(MoS2_lattice, MoS2_crystal,distribution_parameters,prob)
        
        if i%100 == 0:
            MoS2_lattice.plot_lattice(False,dst_data,MoS2_lattice.time[-1],i)
            print (i,MoS2_lattice.time[-1])
            
            #Proportion between adsort 
            #print(i, MoS2_lattice.n_defects,prob[9])
