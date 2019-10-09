# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 15:22:12 2019

@author: guila
"""
import numpy as np
from functions.mass_functions import  A2D, D2R
from shutil import copyfile
import os
from sympy.solvers import solve
from sympy import Symbol

##############################################################################

# This script can be run directly once MDA/MDO has been run.It will run the 
# OpenVSP files automatically as well.

##############################################################################

def geometry_generator():
    print('\n')
    print('calculating the geometry ...')
    print('-- first stage--')
    
###############################################################################
##################          FIRST STAGE   #####################################
###############################################################################
    

    L1 = 9 # initialization of the tanks
    L2 = 21
    
    # get values from OpenMDAO results
    csvData = np.genfromtxt('geom_stage1_file.csv', delimiter=',')
    surfaceData = np.genfromtxt('surface_file.csv', delimiter=',')
    A_e = csvData[1]
    A_t = csvData[2]
    m_add = csvData[3] 
    m_t_lox = csvData[4]
    m_t_h2 = csvData[5]   
    v_lox = csvData[6]
    v_h2 = csvData[7] 
    p_lox = csvData[8]
    p_h2 = csvData[9]
    m_n = csvData[10]
    A_c = csvData[11]
    L_n = csvData[12]
    m_cc = csvData[13]
    l_c = csvData[14]
    h_cone = csvData[15]
    e_n = csvData[16]
    e_cc =  csvData[17]
    A_s = surfaceData[0]
    
    # nozzle geometry
    e_nozzle_1 = e_n
    h_nozzle_1 = L_n
    D_throat_1 = A2D(A_t)
    D_exit_1 = A2D(A_e)
    
    # combustion chamber geometry
    e_cc_1 = e_cc
    h_cone_1 = h_cone
    D_cc_1 = A2D(A_c)
    h_cc_1 =l_c- h_cone_1 
    
    # stage calculation
    geometry_stage_1 = geom_stage(m_add, A_s, L1, L2)
    e_stage_1 = geometry_stage_1[0]
    L_stage_1 = geometry_stage_1[1]
    D_stage_1 = geometry_stage_1[2]
    D_tank_1 = 0.8 * D_stage_1
    A_tank_1 = np.pi * D_tank_1 **2 /4
    
    # LOX tank calculation
    geometry_tank_lox_1 =  geom_tanks(m_t_lox, v_lox, A_tank_1,p_lox)
    L_tank_lox_1 = geometry_tank_lox_1[0]
    e_tank_lox_1 = geometry_tank_lox_1[1]
    D_tank_lox_1 = geometry_tank_lox_1[2]
    L_total_lox_1 = geometry_tank_lox_1[3]
    
    # H2 tank calculation
    geometry_tank_h2_1 =  geom_tanks(m_t_h2, v_h2, A_tank_1,p_h2)
    L_tank_h2_1 = geometry_tank_h2_1[0]
    e_tank_h2_1 = geometry_tank_h2_1[1]
    D_tank_h2_1 = geometry_tank_h2_1[2]
    L_total_h2_1 = geometry_tank_h2_1[3]
        

   
    # set errors between internal wall and tanks        
    err1 = -e_stage_1 + D_stage_1/2 - D_tank_lox_1/2 - e_tank_lox_1
    err2 = -e_stage_1 + D_stage_1/2 - D_tank_h2_1/2 - e_tank_h2_1

    # iteration process to get height and diameter of stage, and diamter of tanks
    while abs(err1) > 0.01  and abs(err2) > 0.01:
        
        geom_stage2_1 = geom_stage(m_add, A_s, L_total_lox_1, L_total_h2_1)
        e_stage_1 = geom_stage2_1[0]
        L_stage_1 = geom_stage2_1[1]
        D_stage_1 = geom_stage2_1[2]
        
        new_Dtank_1 = D_tank_lox_1 +  2 * max(e_tank_h2_1,e_tank_lox_1) - 0 # new tank distance
        new_Atank_1 = np.pi * new_Dtank_1 ** 2 / 4
        
        geom_tank_lox2_1 =  geom_tanks(m_t_lox, v_lox, new_Atank_1,p_lox)
        L_tank_lox_1 = geom_tank_lox2_1[0]
        e_tank_lox_1 = geom_tank_lox2_1[1]
        D_tank_lox_1 = geom_tank_lox2_1[2]
        L_total_lox_1 = geom_tank_lox2_1[3]
        
        geom_tank_h2_1 =  geom_tanks(m_t_h2, v_h2, new_Atank_1,p_h2)
        L_tank_h2_1 = geom_tank_h2_1[0]
        e_tank_h2_1 = geom_tank_h2_1[1]
        D_tank_h2_1 = geom_tank_h2_1[2]
        L_total_h2_1 = geom_tank_h2_1[3]
        
        err1 = -e_stage_1 + D_stage_1/2 - D_tank_lox_1/2 - e_tank_lox_1
        err2 = -e_stage_1 + D_stage_1/2 - D_tank_h2_1/2 - e_tank_h2_1

    e_stage_1 = geom_stage2_1[0]
    L_stage_1 = geom_stage2_1[1]
    D_stage_1 = geom_stage2_1[2]

    


    ## to .as file
    geom_array =  [None] * 19

    geom_array[0] = 'double e_tank_lox_1 = ' + str(e_tank_lox_1) + ';\n'
    geom_array[1] = 'double D_tank_lox_1 = ' + str(D_tank_lox_1) + ';\n'
    geom_array[2] = 'double L_tank_lox_1 = '  + str(L_tank_lox_1) +';\n'
    geom_array[3] = 'double L_total_lox_1 = ' + str(L_total_lox_1)+';\n'
    
    geom_array[4] = 'double e_tank_h2_1 = ' + str(e_tank_h2_1)+';\n'
    geom_array[5] = 'double D_tank_h2_1 = ' + str(D_tank_h2_1)+';\n'
    geom_array[6] = 'double L_tank_h2_1 = ' + str(L_tank_h2_1)+';\n'
    geom_array[7] = 'double L_total_h2_1 = ' + str(L_total_h2_1)+';\n'
    
    geom_array[8] = 'double e_stage_1 = ' + str(e_stage_1)+';\n'
    geom_array[9] = 'double L_stage_1 = ' + str(L_stage_1)+';\n'
    geom_array[10] = 'double D_stage_1 = ' + str(D_stage_1)+';\n'
    
    geom_array[11] = 'double e_nozzle_1 = ' + str(e_nozzle_1)+';\n'
    geom_array[12] = 'double h_nozzle_1 = ' + str(h_nozzle_1)+';\n'
    geom_array[13] = 'double D_throat_1 = ' + str(D_throat_1)+';\n'
    geom_array[14] = 'double D_exit_1 = ' + str(D_exit_1)+';\n'
    
    geom_array[15] = 'double e_cc_1 = ' + str(e_cc_1)+';\n'
    geom_array[16] = 'double h_cone_1 = ' + str(h_cone_1)+';\n'
    geom_array[17] = 'double D_cc_1 = ' + str(D_cc_1)+';\n'
    geom_array[18] = 'double h_cc_1 = ' + str(h_cc_1)+';\n'
    
 
    

    file_name = to_OpenVSP_script(geom_array)
    
    #run the (.as) script created
    myCmd1 = 'vsp -script '+ file_name
    os.system(myCmd1)
    myCmd2 = 'vsp LAST.vsp3'
    os.system(myCmd2)
    return


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



def to_OpenVSP_script(geometrylist = []):
    j = 0
    processing_text = True
    tmptxt = 'LAST_inter.as'
    txt = 'LAST.as'
    with open(tmptxt, 'w') as outfile:
        with open(txt, 'r') as f:
            for line in f:
                if line.startswith('//'):
                    processing_text = False
                    outfile.write(line)
                else: 
                        if line.startswith('    double')and processing_text or line.startswith('double') and processing_text:
                            outfile.write(  '    '+geometrylist[j] )
                            j = j +1
                        else:
                                outfile.write(line)
    output = 'LAST.as'
    copyfile(tmptxt, output)
    os.remove("LAST_inter.as")
    print('\n')
    print('OpenVSP script generated: LAST.as')    
    return output






file_name = geometry_generator()
