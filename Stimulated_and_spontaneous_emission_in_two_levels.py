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
 

N1_0 = 1 #Initial excited state transition
N0_0 = 0 #Initial ground state transition
tau = 1 #Natural decay constant of excited state
T=100000 #Temperature
c=spc.c #Speed of light
lam = 600*10**(-9) #Wavelength of transition
nu= c/lam
x=spc.h*nu/(spc.k*T) #coefficient in Plancks#law
B_T = 2*spc.h*(nu**3)/((np.exp(x)-1)*(c**2))
rho_energy = 4*spc.pi*B_T/c

A_21 = 1
B_21 = A_21/rho_energy
B_12 = B_21   #no degeneracy


def f(N,t):
    N0 = N[0]
    N1 = N[1]
    dN0 = N1/tau+A_21*N1+B_21*N1*rho_energy-B_12*N0*rho_energy
    dN1 = -N1/tau-A_21*N1-B_21*N1*rho_energy+B_12*N0*rho_energy
    return [dN0,dN1]

t0 = 0
tEnd = 1
steps = 100000

fig, ax = plt.subplots(1,1,figsize=(8,4))
t = np.linspace(t0,tEnd,steps)

N = spi.odeint(f,[N0_0,N1_0],t)

ax.plot(t,N[:,1],label="Excited state population")
ax.legend()