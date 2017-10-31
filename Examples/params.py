


#(Perplex) P&T parameter space definitions for phase diagram
#reducing resolution can dramatically reduce run time of Perple_x files
#however, this comes at a cost of accuracy in mineralogy
Pressure_range_mantle_UM    = '3000 1400000'
Temperature_range_mantle_UM = '1400 3000'
resolution_UM               = '60 60' 

Pressure_range_mantle_LM    = '1250000 6500000'
Temperature_range_mantle_LM = '2500 5000'
resolution_LM               = '50 50'


#layers, like concentric shells set here in each region: core, mantle, h20 envelope
num_mantle_layers = 1000
num_core_layers   = 1000
number_h2o_layers = 0

#temperature at surface if no water layer. Essentially temperature below the crust
Mantle_potential_temp = 1700.

#h2o potential Temp, surface temperature if there exists an h2o layer
T_surface_h2o = 300. # K

#initialize planet with these guesses for radial fraction of core and water layer
Core_rad_frac_guess = .54
h20_radfrac_guess   = 0.0

#request exoplex to print extra details on progress and results
verbose = True

#run only perple_x solution files?
perplex_only = False

#multiprocess perple_x?
multi_process = True








