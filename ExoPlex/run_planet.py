# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.

import os
import sys
import minphys
import numpy as np
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))

import functions
import run_perplex
########
#Make files only
only = True

def run_planet_radius(radius_planet, compositional_params, structure_params, layers,filename):

    #find compositional percentages: abun. of each element, core mass frac, core composition and Perplex inputs
    Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = functions.get_percents(*compositional_params)

    #(Perplex)Run fine mesh grid, Upper mantle mineralogy
    Mantle_filename = run_perplex.run_perplex(*[Mantle_wt_per,compositional_params,[structure_params[0],structure_params[1],structure_params[2]],filename,True])

    ##store upper mantle data grids: T, P, rho etc.
    grids_low, names = functions.make_mantle_grid(Mantle_filename,True)
    names.append('Fe')

    #if there is a water mass fraction 0, then append h2o phases to phase list
    if layers[-1] > 0:
        names.append('liq_water')
        names.append('ice_VII')
        names.append('ice_VI')
        names.append('ice_Ih')

    #(Perplex) run the lower mantle grid. Coarse mesh, store data in arrays
    Mantle_filename = run_perplex.run_perplex(*[Mantle_wt_per,compositional_params,[structure_params[3],structure_params[4],structure_params[5]],filename,False])
    grids_high = functions.make_mantle_grid(Mantle_filename,False)[0]

    #only make perplex files?
    if only:
        return

    ###Append low and high res grids. grids are rho, alpha, Cp, T,P from perplex solution
    grids = [grids_low,grids_high]

    #find mass of planet as a function of radius and composition
    Planet = functions.find_Planet_radius(radius_planet, core_mass_frac,structure_params, compositional_params, grids, Core_wt_per, layers)


    Planet['mass'] = minphys.get_mass(Planet,layers)

    Planet['phases'] = functions.get_phases(Planet, grids, layers)

    Planet['phase_names'] = names

    Planet['Vphi'], Planet['Vp'], Planet['Vs'] = functions.get_speeds(Planet, Core_wt_per, grids, layers)

    return Planet
