# coding: utf-8
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.


#**********************************************************************#
'''
This module is intended to run Perple_X through a parameter space to
create corresponding phase diagrams for each composition. Perple_X is
computationally expensive and very time consuming so creating these
mineral data files for a range of compositions before running ExoPlex
is a highly recommended.

To use this module properly, the user must input their compositions with
the frame of mind that they are defining a bulk planet composition. The
inputs are convolved and converted into stoichiometric inputs for perple_x.
There are two ways to do this:

1) The user may input bulk planet composition
    - key parameters:
        define_mantle = False
    -keep in mind all composition inputs are global aso the core and
    mantle are convolved
2) Define mantle composition directly
    - key parameters:
        define_mantle = {'fix_man': True, 'wtCore': .50}
        fFeO = 0.0  *must be zero because we are directly handling mantle composition


Beyond composition, use the default parameters which are located in params.py

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


from params import *


######
'''
Bulk elemental abundances input as lists or arrays (e.g. FeMg = (Fe/Mg)_mol)

'''
######

define_mantle = 0.0 #or dictionary to define mantle comp directly: = {'fix_man': True, 'wtCore': .50}

CaMg = [0.0, 0.06]
AlMg = [0.0, 0.08]
FeMg = [0.5,.7,.8,.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7, 1.8, 1.9, 2.0]
SiMg = [0.5,.7,.8,.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8, 1.9, 2.0]

####
'''
global fraction of Fe moles in the mantle. Essentially a global Fe redox
state parameter.
'''
####

fFeO = 0.0

######
'''
Core composition of Si, O, and S input as core mass fraction. 
These have no impact if you are defining the mantle directly. 
Otherwise, they will impact the stoichiometry of the mantle.
'''
######

wt_frac_Si_core = 0. #0 to 30wt% Si in core
wt_frac_O_core  = 0.
wt_frac_S_core  = 0.



# lists of compositional and structural inputs used to build planet

structure_params = [Pressure_range_mantle_UM, Temperature_range_mantle_UM, resolution_UM,
                    Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                    Core_rad_frac_guess, Mantle_potential_temp, h20_radfrac_guess, T_surface_h2o]

layers = [num_mantle_layers, num_core_layers, number_h2o_layers]

Mass = 1.0
for f in range(len(FeMg)):
    for s in range(len(SiMg)):
        for c in range(len(CaMg)):
            for a in range(len(AlMg)):

                compositional_params = [wt_frac_water, FeMg[f], SiMg[s], CaMg[c], AlMg[a], \
                                        fFeO, wt_frac_Si_core, \
                                        wt_frac_O_core, wt_frac_S_core]
                #run routines to build perplex solution
                Planet = exo.run_perplex_only(compositional_params,structure_params,layers,'planet', define_mantle)







