# -*- coding: utf-8 -*-
"""
Created on Mon May 13 10:36:33 2019

@author: guila
"""
import math as mt
import numpy as np
# Convert Area to Diameter
def A2D(A: float):
    d = np.abs(mt.sqrt(4 * A / np.pi))
    return d

# Convert Diameter to Radius
def D2R(d: float):
    r = d / 2
    return r

# Liquid Oxigen propellant mass calculation
def m_plox_calculation(OF: float, m_p: float):
    m_p_lox = m_p * OF /(OF + 1)
    return m_p_lox
    
# Liquid Oxigen propellant mass calculation
def m_ph2_calculation(OF: float, m_p: float):
    m_p_h2 = m_p  /(OF + 1)
    return m_p_h2
    
# Liquid Oxigen volume calculation
def vlox_calculation(rho_lox: float, OF: float, m_p_lox: float):
    v_lox = 1 / rho_lox * m_p_lox  #different in paper
    return v_lox

# Hydrogen volume calculation
def vh2_calculation(rho_h2: float, OF: float, m_p: float):
    v_h2 = 1 / rho_h2 * m_p  ##different in paper
    return v_h2

# Pressure tank calculation
def pt_calculation(volume: float):
    pt =  (10 ** (-0.1068 * (np.log(volume)-0.2588))) * 1e6
    #if pt < 2e5 or pt > 5e5:
    #    pt = 2e5
    
    return pt
