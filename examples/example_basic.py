
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
    FeMg = 1.
    SiMg = 1.
    wt_frac_Si_core = 0.
    wt_frac_water = 0.
    mol_frac_Fe_mantle = 0
    Pressure_range_mantle = '1000 2500000'
    Temperature_range_mantle =  '1400 4000'

    verbose = True

    #Feel free to change, but these are defaulted right now
    #Do you want to pin Mass and solve Radius?


    #SolveRad = True


    #Earth Mass, Earth Radius, array, array, array,2D grid
    #Mass, Radius, pressure, temperature, density, composition =
    exo.run_planet(FeMg,SiMg,0.011/0.165,0.015/.165,mol_frac_Fe_mantle,wt_frac_Si_core, 0, 0, wt_frac_water, \
               Pressure_range_mantle,Temperature_range_mantle,verbose)