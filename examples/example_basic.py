
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))


import ExoPlex as exo

if __name__ == "__main__":

    #First user must input the ratios
    Radius_planet = 1.

    FeMg = 0.758578
    SiMg = 0.933254

    CaMg = 0.0501187
    AlMg = 0.0794328

    wt_frac_Si_core = 0.
    wt_frac_water = 0.
    mol_frac_Fe_mantle = 0.0
    Pressure_range_mantle = '5000 2000000'
    Temperature_range_mantle =  '1400 4000'

    verbose = True

    #Feel free to change, but these are defaulted right now
    #Do you want to pin Mass and solve Radius?
    resolution = '100 125'

    #FeMg = .149 / .165
    #SiMg = 0.158 / 0.198
    #CaMg = 0.013/0.198
    #AlMg = 0.018/0.198


    wt_frac_O_core = 0.
    wt_frac_S_core = 0.

    num_mantle_layers = 1000
    num_core_layers = 1000
    number_h2o_layers = 0

    Mantle_potential_temp = 1700.
    filename = 'HIP99240_'+str(int(Mantle_potential_temp))+'_rad'+str(int(Radius_planet))


    compositional_params = [wt_frac_water,FeMg,SiMg,CaMg,AlMg,mol_frac_Fe_mantle,wt_frac_Si_core, \
                          wt_frac_O_core,wt_frac_S_core]

    Core_rad_frac_guess = .54
    structure_params =  [Pressure_range_mantle,Temperature_range_mantle,resolution,Core_rad_frac_guess,Mantle_potential_temp]

    layers = [num_mantle_layers,num_core_layers,number_h2o_layers]
    Planet = exo.run_planet_radius(Radius_planet,compositional_params,structure_params,layers,filename)

    print
    print "Mass = ", '%.3f'%(Planet['mass'][-1]/5.97e24), "Earth masses"

    print "CMB Pressure",'%.3f'%(Planet['pressure'][num_core_layers+1]/10000.)
    exo.functions.write(Planet,filename)



    figure = plt.figure(figsize = (12,10))
    # figure.suptitle('Your planet is %.3f Earth Masses with Average Density of %.1f kg/m$^3$' %((Plan.mass/5.97e24), \
    #                (Plan.mass/(4./3.*np.pi*Plan.radial_slices[-1]*Plan.radial_slices[-1]*Plan.radial_slices[-1]))),\
    #                 fontsize=20)

    ax1 = plt.subplot2grid((6, 3), (0, 0), colspan=3, rowspan=3)
    ax2 = plt.subplot2grid((6, 3), (3, 0), colspan=3, rowspan=1)
    ax3 = plt.subplot2grid((6, 3), (4, 0), colspan=3, rowspan=1)
    ax4 = plt.subplot2grid((6, 3), (5, 0), colspan=3, rowspan=1)

    ax1.plot(Planet['radius'] / 1.e3, Planet['density'] / 1.e3, 'k', linewidth=2.)
    ax1.set_ylim(0., (max(Planet['density']) / 1.e3) + 1.)
    ax1.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax1.set_ylabel("Density ( $\cdot 10^3$ kg/m$^3$)")

    # Make a subplot showing the calculated pressure profile
    ax2.plot(Planet['radius'] / 1.e3, Planet['pressure'] / 1.e4, 'b', linewidth=2.)
    ax2.set_ylim(0., (max(Planet['pressure']) / 1e4) + 10.)
    ax2.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax2.set_ylabel("Pressure (GPa)")

    # Make a subplot showing the calculated gravity profile
    ax3.plot(Planet['radius'] / 1.e3, Planet['gravity'], 'r', linewidth=2.)
    ax3.set_ylabel("Gravity (m/s$^2)$")
    ax3.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax3.set_ylim(0., max(Planet['gravity']) + 0.5)

    # Make a subplot showing the calculated temperature profile
    ax4.plot(Planet['radius'] / 1.e3, Planet['temperature'], 'g', linewidth=2.)
    ax4.set_ylabel("Temperature ($K$)")
    ax4.set_xlabel("Radius (km)")
    ax4.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax4.set_ylim(0., max(Planet['temperature']) + 100)

    plt.show()
