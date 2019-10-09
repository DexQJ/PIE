# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 15:22:41 2019

@author: guila
"""
import scipy as sp
import constants as cte
from modele_atmos import compute_atmos
from MODULES.TRAJECTORY.command_law_1 import command_law_1
#from MODULES.TRAJECTORY.command_law_1 import command_law_1
def EOM(t,x, omega, Thrust, w, A,Cd, mu,R,integrate,t_a, t_vertical):                   #Equation of motion
    
    # current state
    r = x[0]                                # radius of trajectory
    l = x[1]                                # longitude of trajectory
    vr = x[2]                               # Radial velocity
    vl = x[3]                               # Tangential velocity
    m = x[4]                                # Mass of system
    
    # Compute current atmosphere
    alt = r-R
    (Temp,Pa,rho,c) = compute_atmos(alt)
    #u = (vr**2+(vl-omega * r )**2)**0.5    #Speed of vehicle with respect to air
    u = (vr**2+vl **2)**0.5                 #Speed of vehicle with respect to air
    Ma = u / c                             # Mach number
    #P = P0/rho0*rho                        # Pressure of air at a given altitude (unused)
    
    
    g_current = mu/((alt+R)**2)
    
    #dynamic pressure
    pdyn = 1 / 2 * rho * u **2
    
    # flux
    flux = pdyn * u
    
    #axial load
    nx = (Thrust(t) - pdyn * A * Cd) /(m * g_current)

    # compute current parameters
    
    D = 1/2 * rho * u  * Cd * A  # Drag coeff 
    
    theta = command_law_1(t_a,t_vertical,t)
    # Differential equations
    dr = vr                                                             
    dl = vl/r
    dvr = vl ** 2 / r + Thrust(t) / m * sp.sin(theta) - g_current - 1 / m * D * vr
    dvl = -vl * vr / r + Thrust(t) / m * sp.cos(theta) - 1 / m * D * (vl)
    dm = -Thrust(t)/w
    

    if integrate == True:
        
        return [dr,dl,dvr,dvl, dm]
    else:
        return[nx,Ma,pdyn,flux,rho,theta]
    
    
    
def nomass(t, x): return x[4] # event for propellant mass totally burned
nomass.terminal = True
nomass.direction = -1


def hit_ground(t, x): return x[0]-cte.RT # event for rocket hitting the ground
hit_ground.terminal = True
hit_ground.direction = -1
