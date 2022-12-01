# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:01:13 2022

@author: ALDANADS
"""
from initialization import *
from KMC import KMC
import shelve

save_data = False
save_var = False


adsortion_rate = [0.00005,0.0001,0.00015,0.0002,0.0004, 0.0006, 0.0008, 0.001, 0.002]
E_nucleation = 1.6
E_propagation = 1.3



parameters = [E_nucleation,E_propagation,adsortion_rate]

for n_sim in np.arange(3,len(adsortion_rate)):

    MoS2_lattice, MoS2_crystal,distribution_parameters,paths,rng = initialization(parameters,n_sim,save_data)
    events = [np.zeros(15),np.zeros(15)]
    
    i = 0
    j = 0
    while MoS2_crystal.coverage < 0.17:
        i += 1
        MoS2_lattice, MoS2_crystal,Mo_adatom,events = KMC(MoS2_lattice, MoS2_crystal,distribution_parameters,events,rng)

        if i%200 == 0:
            j += 1
            MoS2_lattice.plot_lattice(False,paths['data'],MoS2_lattice.time[-1],j,False,MoS2_crystal.cluster_ij)
            print ('Step: ',i, ' Time (s): ',round(MoS2_lattice.time[-1],4),' Coverage (%): ',round(100*MoS2_crystal.coverage,4))

# Variables to save
variables = {'MoS2_lattice' : MoS2_lattice, 
           'MoS2_crystal': MoS2_crystal}

if save_var: save_variables(paths['program'],variables)


