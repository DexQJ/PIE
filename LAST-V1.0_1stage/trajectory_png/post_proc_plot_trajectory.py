# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 17:23:06 2019

@author: guila
"""
import matplotlib.pyplot as plt
import numpy as np
##############################################################################


# Run the script if you either want to plot MDA result OR MDO result
# Please note that you are plotting the latest data saved in trajectory_file.csv. 
# You need to be aware of what you are plotting, otherwise trajectory images 
# will be deleted. 


###############################################################################


def plot_trajectory(MDA):
    if MDA == False:
        my_data = np.genfromtxt('trajectory_file.csv', delimiter=',')
        legend = ['time','alt_ascent (km)','V_ascent (m.s-1)','m_ascent (t)','gamma_ascent (deg)','theta_ascent (deg)',
              'pdyn_ascent (kPa)','lat_ascent','longi_ascent','Mach_ascent','flux_ascent',
              'rho_ascent (kg.m-3)','CX_ascent','thrust_ascent (kN)','Nx_ascent','alpha_ascent']
        for i in range(len(my_data)):
            plt.figure(num=i)
            plt.plot(my_data[0], my_data[i])
            plt.ylabel(legend[i])
            plt.xlabel(legend[0])
            plt.title(legend[i] + '_over_'+legend[0] )
            plt.savefig('MDO/' + legend[i] + '_over_'+ legend[0] + '.png')
    else:
        my_data = np.genfromtxt('trajectory_file.csv', delimiter=',')
        legend = ['time','alt_ascent (km)','V_ascent (m.s-1)','m_ascent (t)','gamma_ascent (deg)','theta_ascent (deg)',
              'pdyn_ascent (kPa)','lat_ascent','longi_ascent','Mach_ascent','flux_ascent',
              'rho_ascent (kg.m-3)','CX_ascent','thrust_ascent (kN)','Nx_ascent','alpha_ascent']
        for i in range(len(my_data)):
            plt.figure(num=i)
            plt.plot(my_data[0], my_data[i])
            plt.ylabel(legend[i])
            plt.xlabel(legend[0])
            plt.title(legend[i] + '_over_'+legend[0] )
            plt.savefig('MDA/' + legend[i] + '_over_'+ legend[0] + '.png')
MDA= True
plot_trajectory(MDA)