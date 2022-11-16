# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:01:13 2022

@author: ALDANADS
"""
from initialization import*
from KMC import*

MoS2_lattice, MoS2_crystal,distribution_parameters = initialization()
prob = np.zeros(len(MoS2_lattice.Act_E)+1)

for i in np.arange(100000):
    MoS2_lattice, MoS2_crystal,Mo_adatom,prob = KMC(MoS2_lattice, MoS2_crystal,distribution_parameters,prob)
    
    if i%100 == 0:
        MoS2_lattice.plot_lattice()
        print (i,MoS2_lattice.time[-1])
        
        #Proportion between adsort 
        #print(i, MoS2_lattice.n_defects,prob[9])
