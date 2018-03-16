#(Perplex) P&T parameter space definitions for phase diagram
#reducing resolution can dramatically reduce run time of Perple_x files
#however, this comes at a cost of accuracy in mineralogy
Pressure_range_mantle_UM    = '5000 1250001'
Temperature_range_mantle_UM = '1500 3200'
resolution_UM               = '30 30'

Pressure_range_mantle_LM    = '1250000 8500000'
Temperature_range_mantle_LM = '2500 5200'
resolution_LM               = '30 30'

#delete current mineral physics file and update on subsequent run?
# check True when you change the values above. Otherwise old params will
# remain
update_minphys = False  

#layers, like concentric shells set here in each region: core, mantle, h20 envelope
num_mantle_layers = 1000
num_core_layers   = 1000
number_h2o_layers = 0

#temperature at surface if no water layer. Essentially temperature below the crust
Mantle_potential_temp = 1800.

#h2o potential Temp, surface temperature if there exists an h2o layer
T_surface_h2o = 300. # K

#water mass fraction
wt_frac_water = 0.0


#initialize planet with these guesses for radial fraction of core and water layer
Core_rad_frac_guess = .54
h20_radfrac_guess   = 0.0

#request exoplex to print extra details on progress and results
verbose = True


#multiprocess perple_x?
multi_process = True


#skip model and run perple_x only?
perplex_only = False#True # False



MEarth = 5.972e24  #kg
REarth = 6.371e6   #meters

#example of how to fix core mass so user can define mantle composition seperately

fix_core = {'fix_man': True, 'wtCore': 0.323}






