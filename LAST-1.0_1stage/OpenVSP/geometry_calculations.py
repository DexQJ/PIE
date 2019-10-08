# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 16:03:46 2019

@author: guila
"""

from sympy.solvers import solve
from sympy import Symbol
import numpy as np
from MODULES.MASS_STRUCTURE.functions.mass_functions import  A2D, D2R


def geom_stage(m_stage: float, A_stage: float, L_lox: float, L_h2, material: str = 'aluminium'):
    if material == 'aluminium':
        rho = 2800                  	        # density of the material (kg/m³)
    D_stage = A2D(A_stage)                      # diameter of stage (m)
    L_stage = (L_lox + L_h2 + 1)                # length of stage + 1 is arbitrary (m)    
    e = Symbol('e')                             # thickness of walls (m)
    
    # solve equation to get thickness ( get postive and real solution)
    sol = solve(-2 * np.pi * L_stage * e ** 2 + (np.pi  * D_stage * L_stage +  np.pi * D_stage ** 2 / 2) * e - m_stage / rho, e )
    if sol[0] > 0 and sol[1] >0:
        e_stage = min(sol)
    else:
        if sol[0] < 0 or sol[1] < 0 or sol[1] == 0 or sol[0] == 0:
            e_stage = max(sol)
        else:
            print('error: e<0 ') 

    # output geometry values            
    return e_stage, L_stage, D_stage


def geom_tanks(m_tank: float, V_tank: float, A_tank: float, P_tank: float, material: str = 'aluminium'):
    if material == 'aluminium':
        rho = 2800                  	                # density of the material (kg/m³)
        sigma = 400 * 1e6                               # maximum stress (Pa)
        ns = 2
    R_tank = D2R(A2D(A_tank))                           # radius of tank (m)
    D_tank = A2D(A_tank)                                # diameter of tank (m)
    L_tank = (V_tank - 4 / 3 *np.pi * (D_tank / 2) ** 3 ) * 4 / (np.pi * D_tank ** 2)   # length of tanks (m)
    if L_tank < 0 :
        L_tank = 0 # if negative value of the length make it 0
    e_tank = P_tank * R_tank *ns / sigma
    L_total_tank = L_tank + D_tank

    # output geometry values
    return L_tank, e_tank, D_tank, L_total_tank




    
        

