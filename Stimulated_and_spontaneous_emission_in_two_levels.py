#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 16:22:12 2018

@author: danielr
"""

#In this simulation, I am trying to simulate how the excited
#state population in a two level system behaves in time,
#where the excited state is fully populated at first and we
#have some kind of stimulated emission. We probably need
#the Einstein equations for that
#the radioation is modelled by Planck's law
#and we are looking at spontanuous decay of the exited state as well
#We are looking at a transition at 600nm

 
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import scipy.constants as spc
 

N1_0 = 1
N0_0 = 0
n_photon = 10
tau_n1 = 1
T=300
c=spc.c
nu= c/(600*10**(-9))
x=spc.k*T/(spc.h*nu)
rho = 8*spc.pi*nu*nu/(spc.c**3)
