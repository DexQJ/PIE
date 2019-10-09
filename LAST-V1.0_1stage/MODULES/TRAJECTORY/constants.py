# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 10:46:18 2019

@author: lbrevaul
"""
import numpy as np
import constants as Cst

#Earth radius (m)
RT = 6378137.0

#Earth mass (kg)
MT =5.9736e24

# density at z=0m (kg/m³)
rho0 = 1.225
#Standard acceleration
g0 = 9.80665

# initial longitude, Kourou (deg)
longitude0 = -52.78

#Pressure on the ground (Pa)
Pa_sol = 101325.0

#Earth rotation time 23h 56 min 4.0905 s
earth_rotation_time = 23*60*60+56*60+4.0905

#Gravity constant
G = 6.67408e-11
mu = 398600.4415e9

#Initial state conditions
r0 = RT+100.
V0 = 10.
gamma0 = np.pi/2.

#Stage separation conditions
Duration_stage_separation = 2.

#Modulation of rocket engine conditions
derating_stage_2 = 0.45
derating_stage_1 = 0.25
coeff_losses = 0.1

#payload mass (kg)
m_pl = 8000

# number of engines
n_engines = 1.