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
import Earth_model as em
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

    

def plot(data1):
    
    fig, ax =  plt.subplots(figsize = (15,10))

    #plotting parameters, can change
    plt.rc('font', family='serif')
    lab_size = 23
    tic_size = 18
    ax.set_xlim(0.0625 , 8)
    ax.set_xlabel(r"Mass (M$_\oplus$)", fontsize = lab_size )
    ax.set_ylabel(r"Radius (R$_\oplus$)", fontsize = lab_size)
    ax.tick_params(direction='in', length=6, labelsize = tic_size)
    ax.grid(color='grey', linestyle='-', alpha = 0.4, linewidth=.7)


    #Plot the PREM and Exoplex model
    #ax.plot(depth_prm, rho_dep, label = 'PREM',  lw = 5, ls = '-.', color = 'black')
    ax.plot(data1[:,0], data1[:,-1], label = 'rock', lw = 4, color = 'blue')
    ax.plot(data1[:,0], data1[:,1], label = '100% Fe', lw = 4, color = 'grey')

    plt.legend(loc = 'lower right', fontsize = tic_size)

    #store the figure somewhere?
    #path_to_figs = '/home/alejandro/Documents/M-R Stuff/ThesisFigs'
    #plt.savefig(path_to_figs+'/Earth_v_PREM_describe.png')

    plt.show()

    
    
    
    
    
dat = ZS_model()
plot(dat)
    
    
    
    
    
    
    
    
    
