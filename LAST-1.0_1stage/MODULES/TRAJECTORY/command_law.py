# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 

@author: aurbano
"""
import numpy as np
import scipy as sp

def command_law_1(dt_po,dthe_po,dt_v,dt_d,fi,t):
#    theta = sp.pi/2 * (-t / dt_po +1)                    # pitch angle
    theta = sp.pi/2 * (-t / dt_po +1)                    # pitch angle
    
    return theta

def command_law_2(dt_po,dthe_po,dt_v,dt_d,fi,t):
    if 0<=t and t<=dt_v:
        theta = np.pi/2
    if dt_v < t:
        theta = np.pi/2 * (1 - (t-dt_v)/dt_po)                    
    
    return theta
 
def command_law_3(dt_po,dthe_po,dt_v,dt_d,fi,t):     #dthe_po in rad
    if 0<=t and t<=dt_v:
        theta = np.pi/2
    if dt_v < t:
#        theta = np.pi/2 *  1 /(t_vertical - t_a ) * (t - t_a)                     
        theta = (np.pi/2  - dthe_po*(t-dt_v)/dt_po)                    

    return theta

def command_law_4(dt_po,dthe_po,dt_v,dt_d,fi,t):     #dthe_po in rad
    if 0<=t and t<=dt_v:
        theta = np.pi/2                              #vertical phase
    elif t>dt_v and t<=dt_po+dt_v:
       theta = (np.pi/2  - dthe_po*(t-dt_v)/dt_po)  #pitch over
    elif t>dt_v+dt_po and t <= dt_v + dt_po + 4*dt_d:
        theta = fi + (np.pi/2 - dthe_po -fi )*sp.exp(-(t-dt_v-dt_po)/dt_d)  # gravity turn 1: return to alfa 0 in dt_d
#        theta = fi - dthe_po*sp.exp(-(t-dt_v-dt_po)/dt_d)  # gravity turn 1: return to alfa 0 in dt_d
    else:
        theta = fi
    return theta
       
def command_law_5(dt_po,dthe_po,dt_v,dt_d,fi,t):     #dthe_po in rad
    if 0<=t and t<=dt_v:
        theta = np.pi/2                              #vertical phase
    elif t>dt_v and t<=dt_po+dt_v:
        theta = (fi  - dthe_po*(t-dt_v)/dt_po)  #pitch over
    elif t>dt_v+dt_po and t <= dt_v + dt_po + 3*dt_d:
        theta = fi - dthe_po*sp.exp(-(t-dt_v-dt_po)/dt_d)  # gravity turn 1: return to alfa 0 in dt_d
    else:
        theta = fi
    return theta
