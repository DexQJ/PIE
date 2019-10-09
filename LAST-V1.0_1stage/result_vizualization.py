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
        legend = ['time','alt_ascent (km)','V_ascent (m.s-1)','m_ascent (t)','gamma_ascent (deg)','theta_ascent (deg)',
              'pdyn_ascent (kPa)','lat_ascent','longi_ascent','Mach_ascent','flux_ascent',
              'rho_ascent (kg.m-3)','CX_ascent','thrust_ascent (kN)','Nx_ascent','alpha_ascent']
        for i in range(len(my_data)):
            plt.figure(num=i)
            plt.plot(my_data[0], my_data[i])
            plt.ylabel(legend[i])
            plt.xlabel(legend[0])
            plt.title(legend[i] + '_over_'+legend[0] )
            plt.savefig('trajectory_png/MDO/' + legend[i] + '_over_'+ legend[0] + '.png')
    else:
        my_data = np.genfromtxt('trajectory_png/trajectory_file.csv', delimiter=',')
        legend = ['time','alt_ascent (km)','V_ascent (m.s-1)','m_ascent (t)','gamma_ascent (deg)','theta_ascent (deg)',
              'pdyn_ascent (kPa)','lat_ascent','longi_ascent','Mach_ascent','flux_ascent',
              'rho_ascent (kg.m-3)','CX_ascent','thrust_ascent (kN)','Nx_ascent','alpha_ascent']
        for i in range(len(my_data)):
            plt.figure(num=i)
            plt.plot(my_data[0], my_data[i])
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
    print('pith over time t_a (s) = ', test['t_a'][0])
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
    print('Thrust1 (kN) =',test['Thrust1'][0] / 1e3)
    print('\n')
    print('final altitude (km) =',test['alt_final'][0])
    print('nx_max ', test['nx_max'][0])
    print('Pdyn_max (kPa) ', test['Pdyn_max'][0])
    print('alpha_max (deg) ', test['alpha_max'][0])
    print('\n')
    

   
    