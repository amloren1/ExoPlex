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

'''
#**********************************************************************#




import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import PREM.prem as p

import pdb
  
############################
'''
Plotting vs prem function
density profile of modeled planet plotted agains the PREM
'''
############################


def pltprem(**kwargs):

    #setup plots
    fig, ax =  plt.subplots(figsize = (15,10))

    #plotting parameters, can change
    plt.rc('font', family='serif', weight = 'medium')
    lab_size = 36
    tic_size = 34
    
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
    ax.set_xlim(0., max(depth_prm))
    rho_dep   = prem_dat.get('rho_depth')
    ax.plot(depth_prm, rho_dep, label = 'PREM',  lw = 5, ls = '-.', color = 'black')
    #in the case of multiple planet dictionaries
    if kwargs.get('planet') is not None:
        Planet = kwargs.get('planet')
        depth_ep = np.empty(len(Planet),dtype = object)
        rho_ep = np.empty(len(Planet),dtype = object)
    
        for i in range(len(Planet)):
            #setup plotting data from ExoPlex model
            #pdb.set_trace()
            depth_ep[i] = (Planet[i]['radius'][-1]-Planet[i]['radius'])/1e3
            rho_ep[i]   = Planet[i]['density']/1e3

            
            #Plot the Exoplex model
            if kwargs.get('label') != None:
                ax.plot(depth_ep[i], rho_ep[i], label = kwargs.get('label')[i] , lw = 4)
            else:
                ax.plot(depth_ep[i], rho_ep[i], label = 'ExoPlex' , lw = 4)
        
    plt.legend(loc = 'lower right', fontsize = 32)

    #store the figure somewhere?
    #path_to_figs = '/home/alejandro/Documents/M-R Stuff/ThesisFigs'
    #plt.savefig(path_to_figs+'/Earth_v_PREM_describe.png')

    plt.show()
    
    #pdb.set_trace()

    #send model data to a file?
    #data = np.array([depth_ep, Planet['density']/1e3])
    #data = data.T
    #np.savetxt('Earth_McD_bulk.dat', data, fmt = '%.4f', delimiter = '\t', header = "depth (km)\tdensity (kg/m^3)")

    return



############################
'''
create data file with all posible columns
'''
############################

def writeall(**kwargs):
    from params import num_core_layers
    
    Planet = kwargs.get('planet')
    n      =  len(Planet)
    
    if kwargs.get('file_names') == None:
        names = np.empty(n, dtype='S20')
        
        for i in range(n):
            names[i] = 'planet0{}.dat'.format(i)
            pdb.set_trace()
    elif len(kwargs.get('file_names')) == n:
        names = kwargs.get('file_names')
    else:
        print '\n\n***ERROR: file_name list does not match number of models***'
        print 'Exiting program'
        sys.exit()
    
    for i in range(len(Planet)):
        
            
        mass = Planet[i]['mass'][num_core_layers:]
        rad = Planet[i]['radius'][num_core_layers:]
        rho = Planet[i]['density'][num_core_layers:]
        P   = Planet[i]['pressure'][num_core_layers:]
        T   = Planet[i]['temperature'][num_core_layers:]
    
        dat_header = '{:25}{:25}{:25}{:25}{:25}'.format('mass', 'radius', 'density', 'pres', 'temp')
        phase_header = ''
        
        n_phase = len(Planet[i]['phase_names'])
        for phases in Planet[i]['phase_names']:
            phase_header += '{:10}'.format(phases)
        #pdb.set_trace()    
        #phase_header = '{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}\
        #    {:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}'.format(*Planet[i]['phase_names'])


        dat = np.transpose([mass, rad, rho, P, T])
        phase = Planet[i]['phases'][num_core_layers:]
        kitchen_sink = np.concatenate([dat,phase],axis=1)


        
        np.savetxt(names[i], kitchen_sink , delimiter = ',' ,header = dat_header+phase_header)
        #np.savetxt('Earth_mantle.dat', np.transpose([mass, rad, rho, P, T]), delimiter = ' , ' ,header = dat_header)
    

############################
'''
write to file without phases
'''
############################

def write(**kwargs):
    from params import num_core_layers
    
    Planet = kwargs.get('planet')
    n      =  len(Planet)
    
    if kwargs.get('file_names') == None:
        names = np.empty(n, dtype='S20')
        
        for i in range(n):
            names[i] = 'planet0{}.dat'.format(i)
            pdb.set_trace()
    elif len(kwargs.get('file_names')) == n:
        names = kwargs.get('file_names')
    else:
        print '\n\n***ERROR: file_name list does not match number of models***'
        print 'Exiting program'
        sys.exit()
    
    for i in range(len(Planet)):
        
            
        mass = Planet[i]['mass'][num_core_layers:]
        rad = Planet[i]['radius'][num_core_layers:]
        rho = Planet[i]['density'][num_core_layers:]
        P   = Planet[i]['pressure'][num_core_layers:]
        T   = Planet[i]['temperature'][num_core_layers:]
    
        dat_header = '{:25}{:25}{:25}{:25}{:25}'.format('mass', 'radius', 'density', 'pres', 'temp')
        phase_header = ''
        
        

        dat = np.transpose([mass, rad, rho, P, T])
        #phase = Planet[i]['phases'][num_core_layers:]
        kitchen_sink = np.concatenate([dat,phase],axis=1)


        #if kwargs.get()
        np.savetxt(names[i], dat , delimiter = ',' ,header = dat_header)
        #np.savetxt('Earth_mantle.dat', np.transpose([mass, rad, rho, P, T]), delimiter = ' , ' ,header = dat_header)
    

############################
'''
plot radial density profile
'''
############################
    
   
def pltrho(**kwargs):
    
    #setup plots
    fig, ax =  plt.subplots(figsize = (15,10))

    #plotting parameters, can change
    plt.rc('font', family='serif', weight = 'medium')
    lab_size = 36
    tic_size = 34
    
    ax.set_ylabel("Density (g/cm$^3$)", fontsize = lab_size )
    ax.set_xlabel("Depth (km)", fontsize = lab_size)
    ax.tick_params(direction='in', length=6, labelsize = tic_size)
    ax.grid(color='grey', linestyle='-', alpha = 0.4, linewidth=.7)

    ax.set_xlim(0., max(depth_prm))
    rho_dep   = prem_dat.get('rho_depth')
    #in the case of multiple planet dictionaries, otherwise
        if kwargs.get('planet') is not None:
        Planet = kwargs.get('planet')
        depth_ep = np.empty(len(Planet),dtype = object)
        rho_ep = np.empty(len(Planet),dtype = object)
    
        for i in range(len(Planet)):
            #setup plotting data from ExoPlex model
            #pdb.set_trace()
            depth_ep[i] = (Planet[i]['radius'][-1]-Planet[i]['radius'])/1e3
            rho_ep[i]   = Planet[i]['density']/1e3

            
            #Plot the Exoplex model
            if kwargs.get('label') != None:
                ax.plot(depth_ep[i], rho_ep[i], label = kwargs.get('label')[i] , lw = 4)
            else:
                ax.plot(depth_ep[i], rho_ep[i], label = 'ExoPlex' , lw = 4)
        
    plt.legend(loc = 'lower right', fontsize = 32)

    #store the figure somewhere?
    #path_to_figs = '/home/alejandro/Documents/M-R Stuff/ThesisFigs'
    #plt.savefig(path_to_figs+'/Earth_v_PREM_describe.png')

    plt.show()
    
    #pdb.set_trace()

    #send model data to a file?
    #data = np.array([depth_ep, Planet['density']/1e3])
    #data = data.T
    #np.savetxt('Earth_McD_bulk.dat', data, fmt = '%.4f', delimiter = '\t', header = "depth (km)\tdensity (kg/m^3)")

    return
    
    
    
    
    
