
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.

import os
import sys
import numpy as np
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))


import ExoPlex as exo

if __name__ == "__main__":

    #First user must input the ratios
    Radius_planet = 1.3

    Star = 'Sun'
    CaMg =0.0616595
    SiMg =0.954993
    AlMg = 0.0851138
    FeMg = 0.812831

    wt_frac_Si_core = 0.
    wt_frac_water = 0.
    mol_frac_Fe_mantle = 0.0

    #(Perplex) P&T parameter space definitions for perplex
    Pressure_range_mantle_UM = '3000 1400000'
    Temperature_range_mantle_UM = '1400 3000'
    resolution_UM = '60 60'

    Pressure_range_mantle_LM = '1250000 6500000'
    Temperature_range_mantle_LM = '2500 5000'
    resolution_LM = '50 50'

    wt_frac_O_core = 0.
    wt_frac_S_core = 0.

    #layers, like concentric shells set here in each region: core, mantle, h20 envelope
    num_mantle_layers = 2000
    num_core_layers = 1000
    number_h2o_layers = 0

    #temperature at surface if no water layer. Essentially temperature below the crust
    Mantle_potential_temp = 1700.

    #h2o potential Temp, surface temperature if there exists an h2o layer
    T_surface_h2o = 300. # K

    #initialize planet with these guesses for radial fraction of core and water layer
    Core_rad_frac_guess = .54
    h20_radfrac_guess = 0.1


    #lists of compositional and structural inputs used to build planet
    compositional_params = [wt_frac_water,FeMg,SiMg,CaMg,AlMg,mol_frac_Fe_mantle,wt_frac_Si_core, \
                          wt_frac_O_core,wt_frac_S_core]


    structure_params =  [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                         Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                         Core_rad_frac_guess,Mantle_potential_temp, h20_radfrac_guess, T_surface_h2o]

    layers = [num_mantle_layers,num_core_layers,number_h2o_layers]
    filename = Star
    Planet = exo.run_planet_radius(Radius_planet,compositional_params,structure_params,layers,filename)

    print
    print "Mass = ", '%.3f'%(Planet['mass'][-1]/5.97e24), "Earth masses"
    print "Core Mass Fraction = ", '%.3f'%(100.*Planet['mass'][num_core_layers]/Planet['mass'][-1])
    print "Core Radius Fraction = ", '%.3f'%(100.*Planet['radius'][num_core_layers]/Planet['radius'][-1])
    print "CMB Pressure = " ,'%.3f' % (Planet['pressure'][num_core_layers]/10000), "GPa"

