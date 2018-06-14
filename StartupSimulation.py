#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 16:06:00 2018

@author: danielr
"""

 #This is the first ever simulation in what I hope
 #will someday become a comprehensive and helpful tool
 #for simulating a lot of Laser physics
 
 #This will be shared on github sometime later, but first I have to build
 #something which is even remotely presentable
 #This first simulation should help me gain experience rapidly
 #It just works by simulating and ploting an exponential decay in a two level
 #system
 
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
 
N0 = 1
tau = 1
t0 = 0
tEnd = 10
steps = 100000
 
 
 
def dN1dt(N1,t):
    return -N1/tau
 
fig, ax = plt.subplots(1,1,figsize=(8,4))
t = np.linspace(t0,tEnd,steps)

N = spi.odeint(dN1dt,N0,t)

ax.plot(t,N,label="Population")
ax.legend()
     
 
 
 