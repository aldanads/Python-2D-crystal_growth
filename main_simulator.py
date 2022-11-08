# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:01:13 2022

@author: ALDANADS
"""
from initialization import*
from KMC import*

MoS2_lattice = initialization()

MoS2_lattice, allowed_events = KMC(MoS2_lattice)

