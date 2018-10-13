#This is a rate-equation simulation of a 4-level Laser system
#The short-lived levels 1 and 3 are ignored, since their population is always nearly zero


from math import *

import scipy.constants as sc
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import rcParams

N0_0, N2_0, phi_0 = 10 ** (20), 0, 0  # Initial conditions: Population in cm^(-3), total number of photons
R1, R2 = 1, 0.95  # Mirror Reflectivities
sigma_abs = 7.7 * 10 ** (-20)  # absorption cross section in cm^2
sigma = 28 * 10 ** (-20)  # Emission cross section in cm^2
T_fluor = 250 * 10 ** (-6)  # Fluorescence lifetime in seconds
P = 0.00001  # Pump power in W
A_pump = 10 ** (-7)  # pump cross section in cm^2
lambda_pump = 0.808 * 10 ** (-4) #Pump wavelength in cm
R_p = sigma_abs * lambda_pump * P / (sc.h * (sc.c * 10 ** 2) * A_pump) #Pump rate
L = 1  # Length of the rod in cm
V_pump = A_pump * L  # Rod volume in cm^3
T_RT = 2 * L / (sc.c * 10 ** 2)  # Round-Trip time in seconds


def four_Level_Laser_model(y, t):
    N0, N2, phi = y[0], y[1], y[2]
    dphi = sc.c * 10 ** 2 * sigma * (phi + 1) * N2 + log(R1 * R2) * phi / T_RT
    dN2 = R_p * N0 - sc.c * 10 ** 2 * sigma * phi * N2 / V_pump - N2 / T_fluor
    dN0 = sc.c * 10 ** 2 * sigma * phi * N2 / V_pump + N2 / T_fluor - R_p * N0
    return [dN0, dN2, dphi]


Time = 10 ** (-4) #Simulation time in seconds
timeUnit = 10**(-6) #time unit in seconds
t = np.linspace(0, Time, 10 ** 7)
y0 = [(N0_0), (N2_0), (phi_0)]

y = odeint(four_Level_Laser_model, y0, t)

rcParams['axes.titlepad'] = 20
ax1 = plt.subplot(211)
plt.tick_params(axis='both', direction='in', labelsize=12)
plt.plot(t/timeUnit, y[0:, 2], label='Number of photons')
plt.title("Rate-equation simulation of a 4-Level Laser system", fontsize=18, fontweight="bold")
plt.ylabel(r"Number of photons [1/$\mathregular{mm^3}$]", fontsize=12)
plt.xticks(np.arange(0, (Time/timeUnit)+1, (Time/timeUnit)/10))
plt.setp(ax1.get_xticklabels(), visible=True)
plt.legend(loc='best')

plt.subplot(212, sharex=ax1)
plt.plot(t/timeUnit, y[0:, 0], label='N0')
plt.plot(t/timeUnit, y[0:, 1], label='N2')
plt.xlabel(r"time [$\mathregular{\mu s}$]", fontsize=12)
plt.ylabel(r"Population [1/$\mathregular{cm^3}$]", fontsize=12)
plt.legend(loc='best')
plt.tick_params(axis='both', direction='in', labelsize=12)
plt.setp(ax1.get_xticklabels(), visible=True)
plt.show()
