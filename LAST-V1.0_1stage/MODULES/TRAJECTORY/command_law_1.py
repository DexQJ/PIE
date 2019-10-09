# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 13:16:58 2019

@author: guila
"""
import numpy as np

def command_law_1(t_a,t_vertical,t):
    
    if 0<=t and t<=t_vertical:
        theta = np.pi/2
    if t_vertical < t:
        theta = np.pi/2 *  1 /(t_vertical - t_a ) * (t - t_a)                     # Thrust angle wrt. inertial frame
    
    return theta
        