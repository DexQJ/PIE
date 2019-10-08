# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 15:34:07 2019

@author: guila
"""
import numpy as np
from MODULES.MASS_STRUCTURE.functions.mass_functions import  A2D, D2R, vlox_calculation, vh2_calculation, pt_calculation, m_plox_calculation, m_ph2_calculation


def mass_additional(m_p: float):
    m_add = m_p * (-2.3 * 1e-7 * m_p + 0.07)
    return m_add

def mass_tank(m_p: float, OF: float, A_e: float, material: str = 'alu'):
    # module parameters 
    if material == 'alu':
        rho = 2800                              # density of the material (kg/m³)
        sigma_max = 400 * 1e6                   # maximum stress (Pa)
        #sigma_el = 300 * 1e6                   # elastic stress (Pa)
        #t_min = 1.5 * 1e-3                     # minimum thickness (mm)
        
    rho_lox = 1141                              # density liquid oxygen (kg/m³)
    rho_h2 = 71                                 # density hydrogen (kg/m³)
    k_t = 2.5                                   # shape coeff
    n_s = 4.0                                   # security factor
    g0 = 9.81
    phi = sigma_max / (rho * k_t * n_s * g0) 
        
    # submodel calculations
    
    # mass of propellant (kg)
    m_plox = m_plox_calculation(OF, m_p)
    m_ph2 = m_ph2_calculation(OF, m_p)
    
    # Volume tanks (m³)
    v_lox = vlox_calculation(rho_lox, OF, m_plox)
    v_h2 = vh2_calculation(rho_h2, OF, m_ph2)
            
    # tank pressures (Pa)
    Pt_lox = pt_calculation(v_lox)
    Pt_h2 = pt_calculation(v_h2)
        

    # output: dry mass of tank
    m_dt = [0] * 6
    m_dt[0] = v_lox * Pt_lox / (g0 * phi)
    m_dt[1] = v_h2 * Pt_h2 / (g0 * phi)
    m_dt[2] = v_lox
    m_dt[3] = v_h2
    m_dt[4] = Pt_lox
    m_dt[5] = Pt_h2
    return m_dt


def mass_turbopump(m_p: float, OF: float, pc: float, m_dot: float):
    # module parameters 
    rho_lox = 1141                              # density liquid oxygen (kg/m³)
    rho_h2 = 71                                 # density hydrogen (kg/m³)
    eta_tp = 0.75                               # turbopump efficiency
    n_lox = 1e4                                 # number of laps per second for LOX
    n_h2 = 3e4                                  # number of laps per second for H2
    A = 1.5                                     # LOX/H2 turbopump coefs (Humble)
    B = 0.6                                     # LOX/H2 turbopump coefs (Humble)
        
    # submodel calculations
    m_p = abs(m_p)
    
    #propellant flow rates (kg/s)
    m_dot_lox = m_dot * OF / (OF + 1)
    m_dot_h2 = m_dot * 1 / (OF + 1)
    
    # propellant mass (kg)
    m_plox = m_plox_calculation(OF, m_p)
    m_ph2 = m_ph2_calculation(OF, m_p)
        
    #volume and pressure of each tank (kg, Pa)
    v_lox = vlox_calculation(rho_lox, OF, m_plox)
    v_h2 = vh2_calculation(rho_h2, OF, m_ph2)
    p_tlox = pt_calculation(v_lox)
    p_th2 = pt_calculation(v_h2)
    
    # power of each turbopump (W)
    w_lox = np.abs(m_dot_lox * (p_tlox - pc) / (rho_lox * eta_tp))
    w_h2 = np.abs(m_dot_h2 * (p_th2 - pc) / (rho_h2 * eta_tp))
    
    # speed of each pump (rad/s)
    omega_lox = 2 * np.pi * n_lox / 60
    omega_h2 = 2 * np.pi * n_h2 / 60
    
    # output: dry mass of turbopump
    m_tp_lox = A * (w_lox / omega_lox) ** B
    m_tp_h2 = A * (w_h2 / omega_h2) ** B
    
    m_dtp = [0] * 2
    m_dtp[0] = m_tp_lox
    m_dtp[1] = m_tp_h2
    
    return m_dtp


def mass_c_chamber(pc: float, A_t: float):
    
    # module parameters
    l_star = 1.02                               # equivalent length for LOx/H2 (Humble)
    n_s = 2.5                                   # security factor
    sigma_s = 310 * 1e6                         # material max stress (Pa)
    rho_s = 8000                                # density of material kg/m³
    theta = 45 * np.pi / 180                    # angle of the cone of the chamber (rad)
  
    # module calculations
    
    # Throat diameter (m)
    d_t = A2D(A_t)
        
    # Chamber Surface (m²)
    A_c = A_t * (8 * (1e2 * d_t) ** (-0.6) + 1.25)
    
    # Chamber diameter (m)
    d_c = A2D(A_c)
    
    # Chamber length (m)
    l_c = l_star * A_t / A_c
    
    # Chamber cone height
    h_cone = (d_c - d_t) / (2 * np.tan(theta))
    
    #chamber thicness
    e_c = pc* d_c * rho_s * n_s / (2 * sigma_s)

    #output: dry mass of the chamber (kg)
    m_cc = (d_c * l_c + (D2R(d_c) ** 2 - D2R(A2D(A_t)) ** 2) / 
            np.sin(theta)) * np.pi * pc* d_c * rho_s * n_s / (2 * sigma_s)
        
    m_dc = np.abs(m_cc)
    return m_dc, l_c, h_cone, A_c,e_c

def mass_nozzle(pc: float, A_e, A_t):
    
    # submodule parameters
    n_s = 2.5                                   # security factor
    sigma_s = 810 * 1e6                         # material max stress(Pa)
    rho_s = 8000                                # density of material kg/m³
    d_e = A2D(A_e)                              # nozzle exit diameter (m)
    d_t = A2D(A_t)                              # throat diameter (m)
    theta_n = 15 * np.pi / 180                  # nozzle exit angle (rad)
    
    # module calculations
    
    # nozzle efficiency (to be included in propulsion module)
    #lambda_n = 1 / 2 * (1 + np.cos(theta_n))
    
    # Chamber Surface (m²)
    A_c = A_t * (8 * (1e2 * d_t) ** (-0.6) + 1.25)
    d_c = A2D(A_c)
    
    # Surface of the nozzle (m²) 
    A_n = np.pi * (D2R(d_e) ** 2 - D2R(d_t) ** 2) / np.sin(theta_n) 
    L_n = (d_e - d_t) / (2 * np.tan(theta_n))
    
    # output: dry mass of nozzle (kg) and length of nozzle (m)
    m_dn = rho_s * A_n * pc * d_c * n_s / (2 * sigma_s) 
    return m_dn, L_n

def mass_press_gaz(OF: float, m_p: float, material: str = 'alu'):
    
    # module parameters        
    if material == 'alu':
        rho = 2800                  	       # density of the material (kg/m³)
        sigma_max = 400 * 1e6                  # maximum stress (Pa)

        
    rho_lox = 1141                              # density liquid oxygen (kg/m³)
    rho_h2 = 71                                 # density hydrogen (kg/m³)
    k_tg = 1.5                                  # shape coefficient
    n_s = 2.0                                   # security factor
    g0 = 9.81                                   # gravity acceleration constant (m/s²)
    
    # module calculations
    
    # gas parameters (He)
    r_g = 2078
    gamma_g = 1.67
    phi_g = sigma_max / (rho * k_tg * n_s * g0) 
    T_gi = 300
    
    # propellant mass (kg)    
    m_plox = m_plox_calculation(OF, m_p)
    m_ph2 = m_ph2_calculation(OF, m_p)
    
    # Volume tanks (m³)
    v_lox = vlox_calculation(rho_lox, OF, m_plox)
    v_h2 = vh2_calculation(rho_h2, OF, m_ph2)
    
    # Pressure of tanks (Pa)
    pt_lox = pt_calculation(v_lox)
    pt_h2 = pt_calculation(v_h2)
    
    # pressurant gas pressure (Pa)
    p_g_lox = 5 * pt_lox 
    p_g_h2 = 5 * pt_h2
    
    # gas mass (kg)
    m_gp_lox = gamma_g * pt_lox * v_lox /(1 - pt_lox / p_g_lox) * (1 / (r_g * T_gi) )        
    m_gp_h2 = gamma_g * pt_h2 * v_h2 /(1 - pt_h2 / p_g_h2) * (1 / (r_g * T_gi) )
        
    # pressurised gas tank mass (kg)
    m_tgp_lox = gamma_g * pt_lox * v_lox / (1 - pt_lox / p_g_lox) * 1 /(g0 * phi_g)
    m_tgp_h2 = gamma_g * pt_h2 * v_h2 / (1 - pt_h2 / p_g_h2) * (1 /(g0 * phi_g))
    
    #output: pressurised gas dry mass
    m_dpg = [0] * 2
    m_dpg[0] = m_gp_lox + m_tgp_lox
    m_dpg[1] = m_gp_h2 + m_tgp_h2
    
    return m_dpg





