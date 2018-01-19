# coding: utf-8
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.


#**********************************************************************#
'''
The goal of this example is to benchmark various models of Earth with
the PREM
The function below, Earth_model, takes in composition and builds an
Earth mass planet. Radius is a function of mass and composition. This can
be inverted as other examples show.
Go from a full description of earth composition and core mass to only
vague knowledge of Earth's bulk elemental abundances indicated by the solar
photosphere.
'''
#**********************************************************************#

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


#prem_dat = 'PREM/PREM.csv'

REarth = 6371 #kilometers

verbose = True


#=======================================================================
'''
this function takes in composition and creates an Earth model

Inputs:
- composition = molar ratios of either bulk planet or just the mantle
                also iron redox state if defining whole planet composition (fFeO)
- corComp     = light element content in the core as mass fraction
- fix_core    = either False or core mass fraction (boolean False or a number between 0-1)
'''
#=======================================================================


def Earth_model(composition, coreComp, fix_core, Mass):
    #import parameters not related to planet composition


    #mass of planet in terms of Earth mass
    #Mass = 1.0 #Earth!

    #Molar ratios for mantle or whole planet
    SiMg       = composition.get('SiMg')
    FeMg       = composition.get('FeMg')
    CaMg       = composition.get('CaMg')
    AlMg       = composition.get('AlMg')
    fFeO_m     = composition.get('fFeO')
    wt_h2o     = composition.get('wt_h2o')


    #composition of core
    wt_frac_Si_core = coreComp.get('wtSi')
    wt_frac_O_core  = coreComp.get('wtO')
    wt_frac_S_core  = coreComp.get('wtS')


    wt_frac_water      = wt_h2o # Earth = 0.0002

    #lists of compositional and structural inputs used to build planet
    compositional_params = [wt_frac_water,FeMg,SiMg,CaMg,AlMg,fFeO_m ,wt_frac_Si_core, \
                          wt_frac_O_core,wt_frac_S_core]

    structure_params =  [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                         Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                         Core_rad_frac_guess,Mantle_potential_temp, h20_radfrac_guess, T_surface_h2o]

    layers = [num_mantle_layers,num_core_layers,number_h2o_layers]

    sol_filename = 'Star_Boy37'


    #Here is where the ExoPlex model is called
    #result is a profile for density, mass, radius,
    #heat capacity, emissivity of heat and mineralogy
    #run_planet_mass(mass_planet, compositional_params, structure_params, layers,filename, truncate_comp)
    
    #Planet = exo.run_planet_mass(Mass, compositional_params,structure_params,layers,sol_filename,  fix_core)

    #run_planet_radius(radius_planet, compositional_params, structure_params, layers,filename, truncate_comp)
    Planet = exo.run_planet_radius(1.0, compositional_params,structure_params,layers,sol_filename, fix_core)
    #print this stuff to make sure you are not going insane in da membrane
    if verbose:

        print 'radius of planet'
        print Planet['radius'][-1]/1000

        print
        print "Mass = ", '%.3f'%(Planet['mass'][-1]/5.97e24), "Earth masses"
        print "Core Mass Fraction = ", '%.3f'%(100.*Planet['mass'][num_core_layers]/Planet['mass'][-1])
        print "Core Radius Fraction = ", '%.3f'%(100.*Planet['radius'][num_core_layers]/Planet['radius'][-1])
        print "CMB Pressure = " ,'%.3f' % (Planet['pressure'][num_core_layers]/10000), "GPa"
        print "Central pressure = {} GPa".format(Planet['pressure'][0]/10000)


    return Planet
    
############################
'''
Plotting function
'''
############################


def plot_vs_PREM(Planet, Planet1):
    
    #import prem data (rad, depth, rho_depth, rho_rad)
    #keys: 'radius', 'depth', 'rho_depth', 'rho_radius' , \
    #     'VPV', 'VPH', 'VSV', 'VSH'}
    prem_dat = p.prem()

    depth_prm = prem_dat.get('depth')
    rho_dep   = prem_dat.get('rho_depth')

    #setup plotting data from ExoPlex model
    depth_ep = (Planet['radius'][-1]-Planet['radius'])/1e3
    rho_ep = Planet['density']/1e3

    depth_ep_Fe = (Planet1['radius'][-1]-Planet1['radius'])/1e3
    rho_ep_Fe = Planet1['density']/1e3


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


    #Plot the PREM and Exoplex model
    ax.plot(depth_prm, rho_dep, label = 'PREM',  lw = 5, ls = '-.', color = 'black')
    ax.plot(depth_ep, rho_ep, label = 'ExoPlex', lw = 4, color = 'magenta')
    ax.plot(depth_ep_Fe, rho_ep_Fe, label = 'ExoPlex Pure Fe core', lw = 4, color = 'green', alpha = 0.5)

    plt.legend(loc = 'lower right', fontsize = 32)

    #store the figure somewhere?
    #path_to_figs = '/home/alejandro/Documents/M-R Stuff/ThesisFigs'
    #plt.savefig(path_to_figs+'/Earth_v_PREM_describe.png')

    plt.show()
    
    pdb.set_trace()

    #send model data to a file?
    data = np.array([depth_ep, Planet['density']/1e3])
    data = data.T
    #np.savetxt('Earth_McD_bulk.dat', data, fmt = '%.4f', delimiter = '\t', header = "depth (km)\tdensity (kg/m^3)")

    return



############################
'''
Print result for modeled planets
'''
############################

def verbo(Planet):
    
        print 'radius of planet in km = {}'.format(Planet['radius'][-1]/1000.)
        print 

        print
        print "Mass = ", '%.3f\n'%(Planet['mass'][-1]/5.97e24), "Earth masses"
        print "Core Mass Fraction = ", '%.3f'%(100.*Planet['mass'][num_core_layers]/Planet['mass'][-1])
        print "Core Radius Fraction = ", '%.3f'%(100.*Planet['radius'][num_core_layers]/Planet['radius'][-1])
        print "CMB Pressure = " ,'%.3f\n' % (Planet['pressure'][num_core_layers]/10000), "GPa"
        print "Central pressure = {} GPa".format(Planet['pressure'][0]/10000)
        
        
        
        
        
        
        
        

#----------------------------------------------------------------

#----------------------------------------------------------------
'''
Suggested composition inputs for the various cases above
'''
#----------------------------------------------------------------

##
#bulk earth ratios from McDonough O3.
# fFeO = 0 and fFeO =  0.1324012
#where fFeO is freaction of total Fe moles in mantle as FeO
##
bulk_Earth_comp = {'FeMg': 2.969696 , 'SiMg': 0.90303030 , \
                 'AlMg': 0.090909090 , 'CaMg': 0.06666 , \
                   'fFeO': 0.0, 'wt_h2o': 0.0}

bulk_inputs = False


bulk_Earth_comp_fFeO = {'FeMg': 2.969696 , 'SiMg': 0.90303030 , \
                 'AlMg': 0.090909090 , 'CaMg': 0.06666 ,  \
                        'fFeO': 0.13240, 'wt_h2o': 0.0}

##
#arth mantle elemental molar ratios from McDonough O3
#fFeO already considered in ratios
# -see man_only parameter below
##
mantle_Earth_comp = {'FeMg': 0.121212121 , 'SiMg': 0.79797979797,  \
              'AlMg': 0.09090909 , 'CaMg': 0.0656565, \
              'fFeO': 0.0,  'wt_h2o': 0.0}


##
#Solar values from Lodders 2003
#use Al/Mg, Ca/Mg composition for Earth mantle but no FeO
##
solar_comp = {'FeMg': 0.813 , 'SiMg': 0.955 , \
                'AlMg': 0.090909090 , 'CaMg': 0.0656565 , \
                'fFeO': 0.0, 'wt_h2o': 0.0}
solar_input = {'fix_man': False, 'wtCore': None}
Fe_only_core = {'wtSi': 0.0, 'wtO': 0.0, 'wtS':0.0}



##
#other composition inputs
# -Earth core composition as wt% of Si, S,O from McDonough 03
# -parameter that tells exoplex to fix core mass and constrain mantle
#  composition seperately
##
light_core_composition = {'wtSi': 0.06, 'wtO': 0.0, 'wtS':0.019}
man_only = {'fix_man': True, 'wtCore': 0.323}

wtCore = .323 #McDonough 03





'''
Run ExoPlex for a variety of cases
'''
import pdb





#########################################################################
#1) Model Earth with full knowledge of composition and core mass fraction

def run():
    plan = Earth_model(mantle_Earth_comp, light_core_composition, man_only, 1.0)

    mass = plan['mass'][num_core_layers:]
    rad = plan['radius'][num_core_layers:]
    rho = plan['density'][num_core_layers:]
    P   = plan['pressure'][num_core_layers:]
    T   = plan['temperature'][num_core_layers:]

    print "***\n\ntest core size and composition of earth analog here\n\n***"
    #pdb.set_trace()
    #dat_header = '{:25}{:25}{:25}{:25}{:25}'.format('mass', 'radius', 'density', 'pres', 'temp')
    #phase_header = '{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}\
     #   {:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}'.format(*plan['phase_names'])


    dat = np.transpose([mass, rad, rho, P, T])
    phase = plan['phases'][num_core_layers:]
    kitchen_sink = np.concatenate([dat,phase],axis=1)

    plan2 = Earth_model(mantle_Earth_comp, Fe_only_core, man_only, 1.0)

    
    #np.savetxt('Earth_nofmt.dat', kitchen_sink , delimiter = ',' ,header = dat_header+phase_header)
    #np.savetxt('Earth_mantle.dat', np.transpose([mass, rad, rho, P, T]), delimiter = ' , ' ,header = dat_header)
    plot_vs_PREM(plan, plan2)

run()

#2) Earth with knowledge of its bulk composition only
#use McD values with and without Fe in mantle
#left light core just because
#Earth_fix_Mcor(bulk_Earth_comp_fFeO, Fe_only_core, bulk_inputs)


#3) Earth mass with solar composition +- errors
#Earth_fix_Mcor(solar_comp, Fe_only_core, solar_input)



##################
