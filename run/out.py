# coding: utf-8
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.


#**********************************************************************#
'''
This is the scripts for all of the native outputs of ExoPlex

contains:
-plot vs prem: density profile of your model vs the prem 
-density profile plot
-mass-radius diagram maker
-table out of all data

TODO:
-plot vs. prem, make sure it works. The idea is to run this fucker from
the scripting file

'''
#**********************************************************************#




import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import PREM.prem as p


  
############################
'''
Plotting function
'''
############################


def plot_vs_PREM(**kwargs):
            
    #setup plots
    fig, ax =  plt.subplots(figsize = (15,10))

    #plotting parameters, can change
    plt.rc('font', family='serif', weight = 'medium')
    lab_size = 36
    tic_size = 34
    ax.set_xlim(0., max(depth_prm))
    ax.set_ylabel("Density (g/cm$^3$)", fontsize = lab_size )
    ax.set_xlabel("Depth (km)", fontsize = lab_size)
    ax.tick_params(direction='in', length=6, labelsize = tic_size)
    ax.grid(color='grey', linestyle='-', alpha = 0.4, linewidth=.7)


    
    #import prem data (rad, depth, rho_depth, rho_rad)
    #keys: 'radius', 'depth', 'rho_depth', 'rho_radius' , \
    #     'VPV', 'VPH', 'VSV', 'VSH'}
    #then plot prem data below
    prem_dat = p.prem()

    depth_prm = prem_dat.get('depth')
    rho_dep   = prem_dat.get('rho_depth')
    ax.plot(depth_prm, rho_dep, label = 'PREM',  lw = 5, ls = '-.', color = 'black')
    #in the case of multiple planet dictionaries
    if kwargs.get('planets') != None:
        Planet = kwargs.get('planets')
        depth_ep = np.empty(len(Planet),dtype = object)
        rho_ep = np.empty(len(Planet),dtype = object)
    
        for i in range(len(Planet)):
        #setup plotting data from ExoPlex model
        depth_ep[i] = (Planet[i]['radius'][-1]-Planet['radius'])/1e3
        rho_ep[i]   = Planet[i]['density']/1e3


    #Plot the PREM and Exoplex model
    
    ax.plot(depth_ep, rho_ep, label = 'ExoPlex', lw = 4, color = 'magenta')
    ax.plot(depth_ep_Fe, rho_ep_Fe, label = 'ExoPlex Pure Fe core', lw = 4, color = 'green', alpha = 0.5)

    plt.legend(loc = 'lower right', fontsize = 32)

    #store the figure somewhere?
    #path_to_figs = '/home/alejandro/Documents/M-R Stuff/ThesisFigs'
    #plt.savefig(path_to_figs+'/Earth_v_PREM_describe.png')

    plt.show()
    
    #pdb.set_trace()

    #send model data to a file?
    data = np.array([depth_ep, Planet['density']/1e3])
    data = data.T
    #np.savetxt('Earth_McD_bulk.dat', data, fmt = '%.4f', delimiter = '\t', header = "depth (km)\tdensity (kg/m^3)")

    return
