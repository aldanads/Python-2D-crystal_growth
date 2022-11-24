# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:01:13 2022

@author: ALDANADS
"""
from initialization import*
from KMC import*
import shutil

#adsortion_rate = [0.0005,0.001,0.002,0.003]
save_data = False

E_nucleation = [1.7, 1.7, 1.7, 1.7, 1.6, 1.6, 1.6, 1.5, 1.5, 1.5]
E_propagation = [1.3, 1.4, 1.5, 1.6, 1.4, 1.5, 1.6, 1.4, 1.5, 1.6]

E_nuc_prop = [E_nucleation,E_propagation]

for n_sim in np.arange(0,len(E_nucleation)):

    MoS2_lattice, MoS2_crystal,distribution_parameters,dst_data = initialization(E_nuc_prop,n_sim,save_data)
    prob = np.zeros(15) # Probability of events
    
    i = 0
    while MoS2_crystal.coverage < 0.25:
        i += 1
        MoS2_lattice, MoS2_crystal,Mo_adatom,prob = KMC(MoS2_lattice, MoS2_crystal,distribution_parameters,prob)

        if i%1 == 0:
            MoS2_lattice.plot_lattice(False,dst_data,MoS2_lattice.time[-1],i,False,MoS2_crystal.cluster_ij)
            print ('Step: ',i, ' Time (s): ',round(MoS2_lattice.time[-1],4),' Coverage (%): ',round(100*MoS2_crystal.coverage,4))
        #break
    #break
            #Proportion between adsort 
            #print(i, MoS2_lattice.n_defects,prob[9])
