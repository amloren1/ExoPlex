#**********************************************************************#
'''
This script is here to make plots comparing ExoPlex with models made by other groups
'''
#*********************************************************************



import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))


import ExoPlex as exo
import PREM.prem as p
from params import *
import pdb
import Earth_models as em
################################################################






####
#   Import the Zeng model
####

def ZS_model():
    
    core_mantle = np.arange(0,23)
    print core_mantle
    f_name = 'Dat/MR_Zeng_short.txt'
    
    
    Fe_rocky = np.genfromtxt(f_name, comments = '#', missing_values=99999, usecols = core_mantle)
    
    
    
    return(Fe_rocky)

    

######################################################
#   Exoplex MvR model for earth mantle and varying core sizes
######################################################
    
    #Earth_model(composition, coreComp, fix_core, Mass)
    
def ExoPlex_model(wtCore):
    
    #Change file name appropriately. Earth core refers to Fe core with 
    # McD 03 light elements amounts
    
    file_name = '{}%_Earth_Core'.format(wtCore*100)
    header = '{:6}{:6}'.format('Mass (ME)', 'radius (RE)')
    
    
    mantle_Earth_comp = {'FeMg': 0.121212121 , 'SiMg': 0.79797979797,  \
              'AlMg': 0.09090909 , 'CaMg': 0.0656565, \
              'fFeO': 0.0,  'wt_h2o': 0.0}

    light_core_composition = {'wtSi': 0.06, 'wtO': 0.0, 'wtS':0.019}
    man_only = {'fix_man': True, 'wtCore': wtCore}

    res  = 50
    dMin = 0.1
    dMax = 5.1
    dMas = (dMax-dMin)/(res-1)
    
    mass_grid   = np.zeros(res)
    radius_grid = np.zeros(res)
    
    for i in range(res):
        mass = dMin+i*dMas
    
        planet = em.Earth_model(mantle_Earth_comp, light_core_composition, man_only, mass)

        mass_grid[i] = planet.get('mass')[-1]/MEarth
        radius_grid[i] = planet.get('radius')[-1]/REarth
    
    
    np.savetxt(file_name, np.transpose([mass_grid,radius_grid]),delimiter = ',',  fmt = '%-10.4f', header = header)
        #pdb.set_trace()
    return(mass_grid, radius_grid)
    

######################################################
#   Plots of M-R
######################################################




def plot(data1):
    
    
    cmf_grid = [0.3,0.5, 0.7]
    
    fig, ax =  plt.subplots(figsize = (15,10))

    #plotting parameters, can change
    plt.rc('font', family='serif')
    lab_size = 23
    tic_size = 18
    ax.set_xlim(0.1 , 5)
    ax.set_xlabel(r"Mass (M$_\oplus$)", fontsize = lab_size )
    ax.set_ylabel(r"Radius (R$_\oplus$)", fontsize = lab_size)
    ax.tick_params(direction='in', length=6, labelsize = tic_size)
    ax.grid(color='grey', linestyle='-', alpha = 0.4, linewidth=.7)


    #Plot the PREM and Exoplex model
    #ax.plot(depth_prm, rho_dep, label = 'PREM',  lw = 5, ls = '-.', color = 'black')
    ax.plot(data1[:,0], data1[:,-1], label = 'rock', lw = 4, color = 'black', ls = ':')
    ax.plot(data1[:,0], data1[:,1], label = '100% Fe', lw = 4, color = 'black', ls = ':')
    ax.plot(data1[:,0], data1[:,11], label = '50% Fe', lw = 4, color = 'black', ls = ':')


    for i in range(len(cmf_grid)):
        #continue
        mas, rad = ExoPlex_model(cmf_grid[i])
        
        ax.plot(mas, rad, label = '{} Core'.format(cmf_grid[i]*100), lw = 4)

    

    plt.legend(loc = 'lower right', fontsize = tic_size)

    #store the figure somewhere?
    #path_to_figs = '/home/alejandro/Documents/M-R Stuff/ThesisFigs'
    #plt.savefig(path_to_figs+'/Earth_v_PREM_describe.png')

    plt.show()


    
dat = ZS_model()
plot(dat)
    
    
    
    
    
    
    
    
    
