# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.

import os
import sys
import minphys
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))
from scipy import interpolate

import functions
import run_perplex
import planet

def run_planet(mass_planet,compositional_params,structure_params,verbose):

    Core_wt_per, Mantle_wt_per, Core_mol_per,core_mass_frac = functions.get_percents(*compositional_params)

    #Mantle_filename,Phases = run_perplex.run_perplex(*[Mantle_wt_per,compositional_params,structure_params])
    Mantle_filename = '0,798_0,903_0,066_0,091_0,0_0,0_MANTLE'
    grids,names = functions.make_mantle_grid(Mantle_filename)

    Planet = planet.initialize(*[mass_planet,core_mass_frac,structure_params,compositional_params])

    Planet = planet.compress(*[Planet,grids,Core_wt_per,structure_params])
    Planet['mass'] = minphys.get_mass(Planet)


    Planet['phases'] = functions.get_phases(Planet,grids,structure_params)

    Planet['phase_names'] = names

    Planet['Vphi'],Planet['Vp'],Planet['Vs'] = functions.get_speeds(Planet,grids,structure_params)


    return Planet
