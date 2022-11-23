# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:01:13 2022

@author: ALDANADS
"""
from initialization import*
from KMC import*
import shutil

#adsortion_rate = [0.0005,0.001,0.002,0.003]
E_nucleation = [1.7, 1.7, 1.7, 1.7, 1.6,1.6,1.6, 1.5, 1.5, 1.5]
E_propagation = [1.3, 1.4, 1.5, 1.6,1.4, 1.5, 1.6, 1.4, 1.5, 1.6]

E_nuc_prop = [E_nucleation,E_propagation]

for n_sim in np.arange(1,2):

    #MoS2_lattice, MoS2_crystal,distribution_parameters,dst_data = initialization(E_nuc_prop,n_sim)
    MoS2_lattice, MoS2_crystal,distribution_parameters = initialization(E_nuc_prop,n_sim)
    prob = np.zeros(len(MoS2_lattice.Act_E)+1)
    
    i = 0
    while MoS2_crystal.coverage < 0.25:
        i += 1
        MoS2_lattice, MoS2_crystal,Mo_adatom,prob = KMC(MoS2_lattice, MoS2_crystal,distribution_parameters,prob)

        if i%100 == 0:
            #MoS2_lattice.plot_lattice(False,dst_data,MoS2_lattice.time[-1],i)
            MoS2_lattice.plot_lattice(False,'',MoS2_lattice.time[-1],i)
            print ('Step: ',i, ' Time (s): ',round(MoS2_lattice.time[-1],4),' Coverage (%): ',round(100*MoS2_crystal.coverage,4))
            
            #Proportion between adsort 
            #print(i, MoS2_lattice.n_defects,prob[9])
