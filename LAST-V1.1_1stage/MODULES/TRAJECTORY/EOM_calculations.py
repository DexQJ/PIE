# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 15:22:41 2019

@author: guila

revised on Mer Aug 28 
@author: aurbano
"""
import scipy as sp
import constants as cte
from modele_atmos import compute_atmos
from MODULES.TRAJECTORY.command_law import command_law_1, command_law_2, command_law_3, command_law_4, command_law_5

def EOM(t,x, omega, ft, w, A, Cd, integrate, dt_po, dt_v, dthe_po,dt_d):                   #Equation of motion
    
    # current state
    r = x[0]                                # radius of trajectory
    l = x[1]                                # longitude of trajectory
    v = x[2]                                # Velocity
    fi = x[3]                               # Flight path angle
    m = x[4]                                # Mass of system
    
    #command law
    theta = command_law_2(dt_po,dthe_po,dt_v,dt_d,fi,t)


    # Compute current atmosphere
    alt = r-cte.RT
    (Temp,Pa,rho,c) = compute_atmos(alt)
    #u = (vr**2+(vl-omega * r )**2)**0.5    #Speed of vehicle with respect to air
 
    # Mach number
    Ma = v / c                              
   
    # gravity acceleration at the specific radius
    g = cte.mu/(r**2)
    
    #dynamic pressure
    pdyn = 1 / 2 * rho * v **2
    
    # thermal flux
    flux = pdyn * v
 
    # drag
    D = pdyn * Cd * A 

    # axial load
    nx = (ft(t) - D * sp.cos(theta - fi)) /(m * cte.g0)


    
    # Differential equations
    dr = v * sp.sin(fi)                                                            
    dl = v * sp.cos(fi)/r
    dv = (ft(t) * sp.cos(theta-fi) - D) / m - g * sp.sin(fi)
    dfi = (v/r - g/v)*sp.cos(fi)+ ft(t)*sp.sin(theta-fi)/ (m*v)
    dm = -ft(t)/w
    

    if integrate == True:
        
        return [dr,dl,dv,dfi,dm]
    else:
        return[nx,Ma,pdyn,flux,rho]
    
    
    
def nomass(t, x): return x[4] # event for propellant mass totally burned
nomass.terminal = True
nomass.direction = -1


def hit_ground(t, x): return x[0]-cte.RT # event for rocket hitting the ground
hit_ground.terminal = True
hit_ground.direction = -1


