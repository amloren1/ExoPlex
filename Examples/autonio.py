import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))


import ExoPlex as exo

run = False
MEarth = 5.972e24 #kg
REarth = 6.371e6   #meters





Star = 'Rando'
CaMg = 0.0 #0.0616595 #
AlMg = 0.0 #0.0851138
FeMg = [0.5,.7,.8,.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7, 1.8, 1.9, 2.0]
SiMg = [0.5,.7,.8,.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8, 1.9, 2.0]


# among light elements, vary only the silicon value
wt_frac_Si_core = 0. #0 to 30wt% Si in core
wt_frac_O_core  = 0.
wt_frac_S_core  = 0.

#h2o content 0 to 1.0
wt_frac_water = 0.0

#how much Fe in the mantle by mole. This is Fe oxidation state
XFeO = 0.0

#(Perplex) P&T parameter space definitions for perplex
#UM-upper mantle, LM-lower mantle
Pressure_range_mantle_UM    = '3000 1400000'
Temperature_range_mantle_UM = '1400 3000'
resolution_UM               = '60 60'

Pressure_range_mantle_LM    = '1250000 6500000'
Temperature_range_mantle_LM = '2500 5000'
resolution_LM               = '50 50'


#layers, like concentric shells set here in each region: core, mantle, h20 envelope
num_mantle_layers = 400
num_core_layers   = 200
number_h2o_layers = 0

#temperature at surface if no water layer. Essentially temperature below the crust
Mantle_potential_temp = 1700.

#h2o potential Temp, surface temperature if there exists an h2o layer
T_surface_h2o = 300. # K

#initialize planet with these guesses for radial fraction of core and water layer
Core_rad_frac_guess = .54
h20_radfrac_guess   = 0.1


structure_params =  [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                     Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                     Core_rad_frac_guess,Mantle_potential_temp, h20_radfrac_guess, T_surface_h2o]









layers = [num_mantle_layers,num_core_layers,0]

Mass = 1.0
for q in range(len(FeMg)):
    for p in range(len(SiMg)):

        compositional_params = [wt_frac_water, FeMg[q], SiMg[p], CaMg, AlMg, XFeO, wt_frac_Si_core, \
                                wt_frac_O_core, wt_frac_S_core]
#run routines to build perplex solution and planet
        Planet = exo.run_planet_mass(Mass, compositional_params,structure_params,layers,'blank', False)







