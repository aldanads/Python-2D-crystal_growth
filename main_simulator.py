# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:01:13 2022

@author: ALDANADS
"""
from initialization import*
from KMC import*

MoS2_lattice = initialization()

for i in np.arange(10):
    MoS2_lattice = KMC(MoS2_lattice)
    MoS2_lattice.plot_lattice()
