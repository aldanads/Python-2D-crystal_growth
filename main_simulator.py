# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:01:13 2022

@author: ALDANADS
"""
from initialization import*
from KMC import*

MoS2_lattice, MoS2_crystal = initialization()

for i in np.arange(30):
    MoS2_lattice, MoS2_crystal = KMC(MoS2_lattice, MoS2_crystal)
    MoS2_lattice.plot_lattice()
