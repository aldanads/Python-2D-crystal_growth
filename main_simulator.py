# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:01:13 2022

@author: ALDANADS
"""
from initialization import*
from KMC import*

MoS2_lattice, MoS2_crystal,pair_atom_defect = initialization()

for i in np.arange(200):
    MoS2_lattice, MoS2_crystal,Mo_adatom = KMC(MoS2_lattice, MoS2_crystal,pair_atom_defect)
    MoS2_lattice.plot_lattice()
