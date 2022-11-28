# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:01:13 2022

@author: ALDANADS
"""
from initialization import*
from KMC import*
import shutil

save_data = True

adsortion_rate = [0.0004,0.0002,0.00001]
E_nucleation = 1.6
E_propagation = 1.3

parameters = [E_nucleation,E_propagation,adsortion_rate]

for n_sim in np.arange(0,len(adsortion_rate)):

    MoS2_lattice, MoS2_crystal,distribution_parameters,dst_data,rng = initialization(parameters,n_sim,save_data)
    events = [np.zeros(15),np.zeros(15)]
    
    i = 0
    j = 0
    while MoS2_crystal.coverage < 0.25:
        i += 1
        MoS2_lattice, MoS2_crystal,Mo_adatom,events = KMC(MoS2_lattice, MoS2_crystal,distribution_parameters,events,rng)

        if i%100 == 0:
            j += 1
            MoS2_lattice.plot_lattice(False,dst_data,MoS2_lattice.time[-1],j,False,MoS2_crystal.cluster_ij)
            print ('Step: ',i, ' Time (s): ',round(MoS2_lattice.time[-1],4),' Coverage (%): ',round(100*MoS2_crystal.coverage,4))
        #break
    #break
            #Proportion between adsort 
            #print(i, MoS2_lattice.n_defects,prob[9])
