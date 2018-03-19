
#**********************************************************************#
'''
This scipt will be for making data files showing a parameter space and
solutions found by ExoPlex
'''
#*********************************************************************



import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))


import ExoPlex as exo
import PREM.prem as p
from params import *
import pdb

################################################################








########################################################################
'''
grids with CMF, FeMg, SiMg input at some mass
'''
########################################################################

def cmf_grid(**kwargs):
    mass = kwargs.get('mass')
    
    #FeMg, SiMg, CMF arrays in form (start,stop, interval)
    FeMg = np.arange(*kwargs.get('femg'))
    SiMg = np.arange(*kwargs.get('simg'))
    CMF  = np.arange(*kwargs.get('cmf'))
    
    if kwargs.get('h2o') == None:
        wt_h2o = [0.0]
    else:
        wt_h2o = np.arange(*kwargs.get('h2o'))
    

    if kwargs.get('filename') == None:
        f_name = '{:2}ME_grid.dat'.format(mass)
    else:
        f_name = '../Solutions/Grids/{}'.format(kwargs.get('filename'))

    
    CaMg       = 0.06
    AlMg       = 0.08
    Fe_ox      = 0.
    wt_Sic, wt_Oc , wt_Sc = 0.,0.,0.
    

    #arrays to store code output for plotting and other nonsense
    #sample exoplex for all these cases mass_grid     = np.zeros(grid_size)
    grid_size = int(len(FeMg)*len(SiMg)*len(CMF)*len(wt_h2o))
    #mass_grid     = np.zeros(grid_size)
    femg_bulk     = np.zeros(grid_size)
    simg_bulk     = np.zeros(grid_size)
    femg_grid     = np.zeros(grid_size)
    simg_grid     = np.zeros(grid_size)
    h2o_grid      = np.zeros(grid_size)
    CMF_grid      = np.zeros(grid_size)
    CRF_grid      = np.zeros(grid_size)
    rad_grid      = np.zeros(grid_size)
    bulk_rho_grid = np.zeros(grid_size)
    n = 0
    for i in range(len(CMF)):
        for j in range(len(FeMg)):
            for k in range(len(SiMg)):
                for h in range(len(wt_h2o)):
                    Fix_core = {'fix_man': True, 'wtCore': CMF[i]}

                    compositional_params = [wt_h2o[h],FeMg[j],SiMg[k],CaMg,AlMg, Fe_ox ,wt_Sic, \
                                  wt_Oc , wt_Sc]

                    structure_params =  [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                                 Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                                 Core_rad_frac_guess,Mantle_potential_temp, h20_radfrac_guess, T_surface_h2o]

                    layers = [num_mantle_layers,num_core_layers,number_h2o_layers]

                    sol_filename = 'NXNN00'

                    Planet = exo.run_planet_mass(mass, compositional_params,structure_params,layers,sol_filename, Fix_core)


                    #mass_grid[n]     = 1.0

                    femg_bulk[n]     = Planet['bulk_ratios'][0]
                    simg_bulk[n]     = Planet['bulk_ratios'][1]
                    femg_grid[n]     = FeMg[j]
                    simg_grid[n]     = SiMg[k]
                    h2o_grid[n]      = wt_h2o[h]
                    CMF_grid[n]      = CMF[i]
                    CRF_grid[n]      = Planet['radius'][num_core_layers]/Planet['radius'][-1]
                    rad_grid[n]      = (Planet['radius'][-1])/REarth #R_Earth units
                    bulk_rho_grid[n] = Planet['mass'][-1]/(4./3 *np.pi*(np.power(Planet['radius'][-1],3)))
                    n+=1

    dat_row_header = '{0:13}{1:13}{2:13}{3:13}{4:13}{5:13}{6:13}{7:13}'.format('CMF','(Fe/Mg)_m','(Si/Mg)_m','Radius',\
                                                '(Fe/Mg)_bulk','(Si/Mg)_bulk','CRF', 'Bulk density')
                                            
    pdb.set_trace()
    print 'header:\n\n'
    print dat_row_header
    np.savetxt(f_name, np.transpose([CMF_grid, femg_grid, simg_grid, rad_grid, femg_bulk,simg_bulk, CRF_grid, bulk_rho_grid]), \
                delimiter = '    ',  fmt = '%-10.4f', header = dat_row_header)

    return

