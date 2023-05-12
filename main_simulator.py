# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:01:13 2022

@author: ALDANADS
"""
from initialization import *
from KMC import KMC
#import shelve

save_data = True
save_var = True


#Temperature = [600, 700, 800, 900, 1000, 1100, 1200]
#etched_adsortion_rate = [0.00005,0.0001,0.00015,0.0002,0.00025,0.0003,0.00035,0.0004,0.00045,0.00050,0.00055,0.0006,0.00065]
non_etched_ad_rate = [0.00015, 0.00021667, 0.00028333, 0.00035, 0.00041667, 0.00048333, 0.00055, 0.00061667, 0.00068333, 0.00075]
etched_adsortion_rate = [0.0002] * len(non_etched_ad_rate)
#non_etched_ad_rate = [0.0002] * 10 
parameters = [non_etched_ad_rate,etched_adsortion_rate]

for n_sim in np.arange(len(parameters[0])):

    MoS2_lattice, MoS2_crystal,distribution_parameters,paths,rng = initialization(parameters,n_sim,save_data)
    events = [0]*15
    
    i = 0
    j = 0
    while MoS2_crystal.coverage < 0.25:
        i += 1
        MoS2_lattice, MoS2_crystal,events = KMC(MoS2_lattice, MoS2_crystal,distribution_parameters,events,rng)

        if i%500 == 0:
            j += 1
            MoS2_lattice.plot_lattice(False,paths['data'],MoS2_lattice.time[-1],j,False,MoS2_crystal.cluster_ij,distribution_parameters[5])
            print ('Step: ',i, ' Time (s): ',round(MoS2_lattice.time[-1],4),' Coverage (%): ',round(100*MoS2_crystal.coverage,4))
            print (sum(np.array(MoS2_lattice.v_adatom_flux)[0])/MoS2_lattice.time[-1]/MoS2_crystal.Mo_sites)
    # Variables to save
    variables = {'MoS2_lattice' : MoS2_lattice, 
               'MoS2_crystal': MoS2_crystal}
    
    if save_var: save_variables(paths['program'],variables)


