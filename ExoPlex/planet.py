import os
import sys
import numpy as np
Earth_radius = 6.371e6
import minphys


def initialize(*args):
    mass_planet = args[0]
    core_mass_frac = args[1]
    structural_params = args[2]
    compositional_params = args[3]
    phases = ['a']

    wt_frac_water, FeMg, SiMg, CaMg, AlMg, mol_frac_Fe_mantle, wt_frac_Si_core, \
     wt_frac_O_core, wt_frac_S_core = compositional_params

    Mantle_potential_temp,num_mantle_layers,num_core_layers,number_h2o_layers = structural_params[3:]


    #def Setup(Mp, masMan, masCor, masH2O, nr, tCrstMan, dtManCor):
    # global Th,Tss
    #    global rad, rho, rho0, vol, dmass, mass, delta,\
    #        T, nic, noc, grav, P, rhon, alpha, cp,phase

    # use inputs to find actual mass of mantle and core

    # if there is a water layer, the imput temperature is lowered because that temperature is for the crustal layer
    # also 50 shells are used for the water layer hence the nh20 vaiable

    if wt_frac_water > 0:

        Surface_temp = 300

    else:
        Surface_temp = Mantle_potential_temp


    if wt_frac_water == 0. and number_h2o_layers > 0:
       print "You have layers of water but no water!"
       number_h2o_layers = 0

    num_layers = num_core_layers+num_mantle_layers + number_h2o_layers  # add 50 shells if there is an h2O layer
    num_phases = len(phases) + 3 #includes ice/water layers
    # arrays to be used

    # these grids are initialized in this function and are then a passed around through the other routines and edited each iteration
    radius_layers = np.zeros(num_layers + 1)
    rho_layers = np.zeros(num_layers + 1)
    rho0_layers = np.zeros(num_layers + 1)
    volume_layers = np.zeros(num_layers + 1)
    mass_layers = np.zeros(num_layers + 1)
    cumulative_mass = np.zeros(num_layers + 1)
    Temperature_layers = np.zeros(num_layers + 1)
    ########################################################################################################## Move to compression function and set as outputs
    # used for compressiofh2on funciton
    gravity = np.zeros(num_layers + 1)
    Pressure_layers = np.zeros(num_layers + 1)

    rhon = np.zeros(num_layers + 1)
    alpha = np.zeros(num_layers + 1)
    cp = np.zeros(num_layers + 1)
    phase = np.zeros((num_layers + 1, num_phases))  # 15 mineral phases + 2ice + liquid water #phasechange
    ############################################################################################################################################

    # Recall that masCor and masMan are decimal percentages
    #mass_water = mass_planet * wt_frac_water
    #Mass_rock_core = mass_planet * (1. - wt_frac_water)  # mass of the planet without h2o layer

    #mass_core = mass_planet * mass_planet  # mass of inner core is half the mass of the core, only one homogenoug core is used now
    #mass_mantle = mass_planet - mass_core

    # setup the initial solution################
    # use data @ T = 2500, P=10,000 fo


    #         @ T = TManCore, P = 10,000 for core
    #
    # initialize plant of pure iron core, pv mantle and LIQUID H2O layer on top



    # core
    # may not have to do this, perplex and other code will take cae of differentiation

    planet_radius_guess = mass_planet*Earth_radius
    core_thickness_guess = (core_mass_frac)*planet_radius_guess
    mantle_thickness_guess = planet_radius_guess - (wt_frac_water*planet_radius_guess) - core_thickness_guess

    water_thickness_guess = planet_radius_guess - mantle_thickness_guess - core_thickness_guess


    for i in range(num_layers+1):
        if i <num_core_layers:
            radius_layers[i]=((float(i)/num_core_layers)*core_thickness_guess)
            rho_layers[i]=8280.0
            rho0_layers[i]=8280.0
            volume_layers[i]=(4.0 / 3.) * np.pi * (pow(radius_layers[i],3.)- pow(radius_layers[i - 1],3.))
            mass_layers[i]=rho_layers[i] * volume_layers[i]
            cumulative_mass[i]=cumulative_mass[i - 1] + mass_layers[i]
            Temperature_layers[i] = 2000.+ Mantle_potential_temp + 0.3 * (core_thickness_guess+mantle_thickness_guess-\
                                                                          - radius_layers[i])/1000.

        elif i <= (num_core_layers+num_mantle_layers):
            radius_layers[i]=(core_thickness_guess+((float(i-num_core_layers)/num_mantle_layers)*mantle_thickness_guess))
            rho_layers[i]=4100.
            rho0_layers[i]=4100.0
            volume_layers[i]= (4.0 / 3.) * np.pi * (pow(radius_layers[i],3.)- pow(radius_layers[i - 1],3.))
            mass_layers[i]= rho_layers[i] * volume_layers[i]
            cumulative_mass[i]= cumulative_mass[i - 1] + mass_layers[i]
            Temperature_layers[i] = Mantle_potential_temp + 0.5 * (core_thickness_guess+mantle_thickness_guess\
                                                                   - radius_layers[i])/1000.


        else:
            print "in here"
            radius_layers[i]=core_thickness_guess+mantle_thickness_guess+\
                             ((float(i-num_core_layers-num_mantle_layers)/number_h2o_layers)*water_thickness_guess)
            rho_layers[i]=1000.
            rho0_layers[i]=1000.
            volume_layers[i]=(4.0 / 3.) * np.pi * (pow(radius_layers[i],3.)- pow(radius_layers[i - 1],3.))
            mass_layers[i]=rho_layers[i] * volume_layers[i]
            cumulative_mass[i]=cumulative_mass[i - 1] + mass_layers[i]
            #Temperature_layers[i] = Surface_temp +(float(i-num_core_layers-num_mantle_layers)/number_h2o_layers)*water_thickness_guess)
            Temperature_layers[i] = 300.

    for i in range(num_layers+1):
        if i >= num_core_layers+num_mantle_layers:
            Pressure_layers[i] = 1
        elif i > num_core_layers:
            Pressure_layers[i] = 10000*(num_core_layers+num_mantle_layers-i)
        else:
            Pressure_layers[i] = 1000000*(num_core_layers-i)


    #initial temperature guess of 0.5 K per km

    return ([radius_layers, rho_layers, rho0_layers, volume_layers,\
             mass_layers, cumulative_mass, Temperature_layers,\
             gravity, Pressure_layers, rhon, alpha, cp, phase])

def compress(*args):
    Planet = args[0]
    Mass_planet = args[1]
    Mantle_filename = args[2]
    compositional_params = args[3]
    structural_params= args[4]


    n_iterations = 1
    max_iterations = 1

    while n_iterations <= max_iterations:
        rho_layers = minphys.get_rho(Planet,Mantle_filename,structural_params)
        n_iterations+=1
