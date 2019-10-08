# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 10:46:18 2019

@author: lbrevaul
"""

import numpy as np
from matplotlib import pyplot as plt


def plot_trajectory(MDA):
    if MDA == False:
        my_data = np.genfromtxt('trajectory_png/trajectory_file.csv', delimiter=',')
        size= len(my_data[1,:])    
        legend = ['time[s]','alt_ascent[km]','V_ascent[m.s-1]','m_ascent[t]','l_ascent[deg]','fi_ascent[deg]','theta_ascent[deg]',
              'alpha_ascent[deg]','thrust_ascent[kN]','pdyn_ascent[kPa]','Mach_ascent','flux_ascent[kW.m-2]',
              'rho_ascent[kg.m-3]','Nx_ascent','Cx_ascent']
        for i in range(size):
            plt.figure(num=i)
            plt.plot(my_data[:,0], my_data[:,i])
            plt.ylabel(legend[i])
            plt.xlabel(legend[0])
            plt.title(legend[i] + '_over_'+legend[0] )
            plt.savefig('trajectory_png/MDO/' + legend[i] + '_over_'+ legend[0] + '.png')
    else:
        my_data = np.genfromtxt('trajectory_png/trajectory_file.csv', delimiter=',')
        size= len(my_data[1,:])        
        legend = ['time[s]','alt_ascent[km]','V_ascent[m.s-1]','m_ascent[t]','l_ascent[deg]','fi_ascent[deg]','theta_ascent[deg]',
              'alpha_ascent[deg]','thrust_ascent[kN]','pdyn_ascent[kPa]','Mach_ascent','flux_ascent[kW.m-2]',
              'rho_ascent[kg.m-3]','Nx_ascent','Cx_ascent']
        for i in range(size):
            plt.figure(num=i)
            plt.plot(my_data[:,0], my_data[:,i])
            plt.ylabel(legend[i])
            plt.xlabel(legend[0])
            plt.title(legend[i] + '_over_'+legend[0] )
            plt.savefig('trajectory_png/MDA/' + legend[i] + '_over_'+ legend[0] + '.png')
    
    
def view_values(P_out):
    test = P_out
    print('\n')
    print('GLOW (t) = ', test['GLOW'])
    print('f1 = ',test['f1'][0])
    print('TW1 = ',test['TW1'][0])
    print('C1 (A_stage / A_exit) = ',test['C1'][0])
    print('pc1 (Pa) = ', test['pc1'][0])
    print('\n')
    print('vertical phase duration dt_v (s) = ', test['dt_v'][0])
    print('pith over duration dt_po (s) = ', test['dt_po'][0])
    print('pith over angle dthe_po (s) = ', test['dthe_po'][0])
    print('return to alf 0 duration (s) = ', test['dt_d'][0])
    print('m_p1 (kg) =',test['m_p_1'][0])
    print('flow rate1 (kg/s) =',test['m_dot1'][0])
    print('\n')
    print('m_dry1 (kg) =',test['m_d_1'][0])
    print('\n')
    print('A_stage1 (m²) =', test['A_s_1'][0])
    print('A_exit1 (m²) = ',test['A_e_1'][0])
    print('A_throat1 (m²) = ',test['A_t_1'][0])
    #print('A_c= ',test['A_c'])
    print('\n')
    print('Isp1 (s) =',test['Isp1'][0])
    print('c*1 (m/s) =',test['c_star1'][0])
    print('epsilon1 =',test['epsilon1'][0])
    print('Thrust1 (kN) =',test['ft1'][0] / 1e3)
    print('\n')
    print('final altitude (km) =',test['alt_f'][0])
    print('nx_max ', test['nx_max'][0])
    print('Pdyn_max (kPa) ', test['Pdyn_max'][0])
    print('alpha_max (deg) ', test['alf_max'][0])
    print('\n')
    

   
    