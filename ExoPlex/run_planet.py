# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.

import os
import sys

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))

import functions as func
import run_perplex

def run_planet(FeMg,SiMg,CaMg,AlMg,mol_frac_Fe_mantle,wt_frac_Si_core, wt_frac_O_core, wt_frac_S_core, wt_frac_water,verbose):

    Core_wt_per, Mantle_wt_per, Core_mol_per = func.get_percents(FeMg,SiMg,0.,0.,0.,0.,0.,0.)

    run_perplex.run_perplex(Mantle_wt_per,SiMg,FeMg,CaMg,AlMg,mol_frac_Fe_mantle,wt_frac_Si_core,wt_frac_O_core,wt_frac_S_core,verbose)




    return 0
