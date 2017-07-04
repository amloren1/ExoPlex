
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
    #Radius_planet = 1.

    wt_frac_Si_core = 0.
    wt_frac_water = 0.
    mol_frac_Fe_mantle = 0.0

    Pressure_range_mantle_UM = '3000 1400000'
    Temperature_range_mantle_UM = '1400 3500'
    resolution_UM = '80 140' #11200

    Pressure_range_mantle_LM = '1250000 6200000'
    Temperature_range_mantle_LM = '2200 5000'
    resolution_LM = '80 80' #6400

    wt_frac_O_core = 0.
    wt_frac_S_core = 0.

    num_mantle_layers = 1500
    num_core_layers = 1000
    number_h2o_layers = 0

    Core_rad_frac_guess = .54

    layers = [num_mantle_layers,num_core_layers,number_h2o_layers]

    Radius_planet = [0.5,.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5]

    Mantle_potential_temp = 1700.

    Mass = []
    CRF = []
    CMF = []

    Star = 'Sun'
    CaMg = 0.0616595
    SiMg = 0.954993
    AlMg = 0.0851138
    FeMg = 0.812831

    compositional_params = [wt_frac_water,FeMg,SiMg,CaMg,AlMg,mol_frac_Fe_mantle,wt_frac_Si_core, \
                          wt_frac_O_core,wt_frac_S_core]

    structure_params = [Pressure_range_mantle_UM, Temperature_range_mantle_UM, resolution_UM,
                        Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                        Core_rad_frac_guess, Mantle_potential_temp]
    for i in Radius_planet:
        filename = Star

        print "Radius", '%.3f' % i

        Planet = exo.run_planet_radius(i,compositional_params,structure_params,layers,filename)
        print
        print "Mass = ", '%.3f' % (Planet['mass'][-1] / 5.97e24), "Earth masses"
        print "Core Mass Fraction = ", '%.3f' % (100. * Planet['mass'][num_core_layers] / Planet['mass'][-1])
        print "Core Radius Fraction = ", '%.3f' % (100. * Planet['radius'][num_core_layers] / Planet['radius'][-1])
        print "CMB Pressure = ", '%.3f' % (Planet['pressure'][num_core_layers] / 10000), "GPa"

        Mass.append(Planet['mass'][-1]/5.97e24)
        CRF.append(100.*Planet['radius'][num_core_layers]/Planet['radius'][-1])
        CMF.append(100.* Planet['mass'][num_core_layers] / Planet['mass'][-1])

        filename = Star+'_'+ str(int(Mantle_potential_temp)) + '_rad_' + str(i)
        exo.functions.write(Planet,filename)

    output_file = Star+"_"+str(int(Mantle_potential_temp))+'_massCRF.dat'
    output = [[Radius_planet[i], Mass[i],CRF[i],CMF[i]] for i in range(len(Mass))]

    np.savetxt(output_file,output , '%.5f', "\t", newline='\n',
               header='Radius\tMass\tCRF\tCMF', footer='', comments='# ')

    #exo.functions.write(Planet,filename)



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

    #plt.show()
