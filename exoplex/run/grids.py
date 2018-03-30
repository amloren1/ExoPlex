
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


import main as exo
import PREM.prem as p
from params import *
import pdb

################################################################







def grid(**kwargs):
    
    ########################################################################
    '''
    grids with range of bulk FeMg, SiMg input at some mass
    '''
    ########################################################################

    mass = kwargs.get('mass')
    
    #FeMg, SiMg, CMF arrays in form (start,stop, interval)
    FeMg = np.arange(*kwargs.get('femg'))
    SiMg = np.arange(*kwargs.get('simg'))
    
    if kwargs.get('h2o') == None or len(kwargs.get('h2o')) == 0:
        wt_h2o = [0.0]
    else:
        wt_h2o = np.arange(*kwargs.get('h2o'))
    

    if kwargs.get('filename') == None:
        f_name = 'Solutions/Grids/{:2}ME_bulk.dat'.format(mass)
    else:
        f_name = 'Solutions/Grids/{}'.format(kwargs.get('filename'))

    
    CaMg       = 0.06
    AlMg       = 0.08
    Fe_ox      = 0.
    wt_Sic, wt_Oc , wt_Sc = 0.,0.,0.
    

    #arrays to store code output for plotting and other nonsense
    #sample exoplex for all these cases mass_grid     = np.zeros(grid_size)
    grid_size = int(len(FeMg)*len(SiMg)*len(wt_h2o))
    #mass_grid     = np.zeros(grid_size)
    femg_man      = np.zeros(grid_size)
    simg_man      = np.zeros(grid_size)
    femg_grid     = np.zeros(grid_size)
    simg_grid     = np.zeros(grid_size)
    h2o_grid      = np.zeros(grid_size)
    CMF_grid      = np.zeros(grid_size)
    CRF_grid      = np.zeros(grid_size)
    rad_grid      = np.zeros(grid_size)
    bulk_rho_grid = np.zeros(grid_size)
    n = 0

    for j in range(len(FeMg)):
        for k in range(len(SiMg)):
            for h in range(len(wt_h2o)):
                Fix_core = {'fix_man': False, 'wtCore': 99999}

                compositional_params = [wt_h2o[h],FeMg[j],SiMg[k],CaMg,AlMg, Fe_ox ,wt_Sic, \
                              wt_Oc , wt_Sc]

                structure_params =  [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                             Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                             Core_rad_frac_guess,Mantle_potential_temp, h20_radfrac_guess, T_surface_h2o]

                layers = [num_mantle_layers,num_core_layers,number_h2o_layers]

                sol_filename = 'NXNN00'

                Planet = exo.run_planet_mass(mass, compositional_params,structure_params,layers,sol_filename, Fix_core)
                

                #mass_grid[n]     = 1.0

                femg_man[n]      = Planet['mantle_ratios'][0]
                simg_man[n]      = Planet['mantle_ratios'][1]
                femg_grid[n]     = FeMg[j]
                simg_grid[n]     = SiMg[k]
                h2o_grid[n]      = wt_h2o[h]
                CMF_grid[n]      = Planet['mass'][num_core_layers]/Planet['mass'][-1]
                CRF_grid[n]      = Planet['radius'][num_core_layers]/Planet['radius'][-1]
                rad_grid[n]      = (Planet['radius'][-1])/REarth #R_Earth units
                bulk_rho_grid[n] = Planet['mass'][-1]/(4./3 *np.pi*(np.power(Planet['radius'][-1],3)))
                n+=1

    dat_row_header = '{0:13}{1:13}{2:13}{3:13}{4:13}{5:13}{6:13}{7:13}'.format('CMF','(Fe/Mg)_bulk','(Si/Mg)_bulk','Radius',\
                                                '(Fe/Mg)_m','(Si/Mg)_m','CRF', 'Bulk density (g/cc)')
                                            
    np.savetxt(f_name, np.transpose([CMF_grid, femg_grid, simg_grid, rad_grid, femg_man,simg_man, CRF_grid, bulk_rho_grid]), \
                delimiter = '    ',  fmt = '%-10.4f', header = dat_row_header)

    print '\n\nGrid stored in {}'.format(f_name)
    print '\nHeader:\n{}'.format(dat_row_header) 
    return





def grid_cmf(**kwargs):
    
    
    ########################################################################
    '''
    grids with CMF, FeMg, SiMg input at some mass
    '''
    ########################################################################

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
        f_name = 'Solutions/Grids/{:2}ME_cmf.dat'.format(mass)
    else:
        f_name = 'Solutions/Grids/{}'.format(kwargs.get('filename'))

    
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
                                                '(Fe/Mg)_bulk','(Si/Mg)_bulk','CRF', 'Bulk density (g/cc)')
                                            
    np.savetxt(f_name, np.transpose([CMF_grid, femg_grid, simg_grid, rad_grid, femg_bulk,simg_bulk, CRF_grid, bulk_rho_grid]), \
                delimiter = '    ',  fmt = '%-10.4f', header = dat_row_header)

    print '\n\nGrid stored in {}'.format(f_name)
    print '\nHeader:\n{}'.format(dat_row_header) 
    return



########################################################################
'''
grids for a range of masses and composition. Compositions are given for BULK PLANET
and core mass fraction is CALCULATED
Each composition is one grid
'''
########################################################################

def mvr_grid(**kwargs):
    
    
    
    #FeMg, SiMg, CMF arrays in form (start,stop, interval)
    mass = np.arange(*kwargs.get('mass'))
    FeMg = np.arange(*kwargs.get('femg'))
    SiMg = np.arange(*kwargs.get('simg'))
    CMF  = np.arange(*kwargs.get('cmf'))
    
    if kwargs.get('h2o') == None:
        wt_h2o = [0.0]
    else:
        wt_h2o = np.arange(*kwargs.get('h2o'))
    
    #filenames come in as lists or arrays

    if kwargs.get('filename') == None:
        default_filename = True
    else:
        default_filename = False
    
   
    
    CaMg       = 0.06
    AlMg       = 0.08
    Fe_ox      = 0.
    wt_Sic, wt_Oc , wt_Sc = 0.,0.,0.
    

    ##arrays to store code output for plotting and file output
    ##sample exoplex for all these cases mass_grid     = np.zeros(grid_size)
    
    #number of rows per file
    n_rows = int(len(mass))
    
    #number of files (each composition gets one)
    n_files = int(len(FeMg)*len(SiMg)*len(CMF)*len(wt_h2o))
    
    rad            = np.zeros((n_files,n_rows))
    rho_bulk       = np.zeros(n_rows)
    CRF_grid       = np.zeros(n_rows)
    CMF_grid       = np.zeros(n_rows)
    
    R     = np.empty(n_files, dtype = object)
    RR     = np.empty(n_files, dtype = object)
    M     = np.empty(n_files, dtype = object)
    R_test = np.empty(0, dtype = object)
    label = np.empty(n_files, dtype = 'S30')
    nn = 0
    
        
    for j in range(len(FeMg)):
        for k in range(len(SiMg)):
            for h in range(len(wt_h2o)):
                
                #set file names for this composition
                if default_filename:
                    f_temp = 'SiMg_FeMg_h2o_bulk{:3}_{:3}_{:3}'.format(SiMg[k], FeMg[j], wt_h2o[h])
                    f_name = 'Solutions/Grids/'+f_temp.replace('.', ',')+'.dat'
                else:
                    f_name = kwargs.get('filename') 
              
                #initialize exoplex inputs for this composition
                Fix_core = {'fix_man': False, 'wtCore': 99999}

                compositional_params = [wt_h2o[h],FeMg[j],SiMg[k],CaMg,AlMg, Fe_ox ,wt_Sic, \
                              wt_Oc , wt_Sc]

                structure_params =  [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                             Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                             Core_rad_frac_guess,Mantle_potential_temp, h20_radfrac_guess, T_surface_h2o]

                layers = [num_mantle_layers,num_core_layers,number_h2o_layers]

                sol_filename = 'NXNN00'

                #run exoplex for range of masses where composition is fixed above
                for m in range(n_rows):

                    #exoplex
                    Planet = exo.run_planet_mass(mass[m], compositional_params,structure_params,layers,sol_filename, Fix_core)
                    
                    rad[nn,m]      = (Planet['radius'][-1])/REarth #R_Earth units
                    rho_bulk[m] = Planet['mass'][-1]/(4./3 *np.pi*(np.power(Planet['radius'][-1],3))) #g/cc
                    CRF_grid[m] = Planet['radius'][num_core_layers]/Planet['radius'][-1]
                    CMF_grid[m] = Planet['mass'][num_core_layers]/Planet['mass'][-1]
                
                #store mass and radius arrays in an object array for plotting
                M[nn]  = mass
                R[nn]  = rad[nn,:]
                #RR[nn] = rad
                label[nn] = 'SiMg={:1},FeMg={:1}'.format(SiMg[k], FeMg[j])
                
                #pdb.set_trace()
                #bulk composition for file header
                femg_bulk     = Planet['bulk_ratios'][0]
                simg_bulk     = Planet['bulk_ratios'][1]
                    
                header = 'Fe/Mg_bulk = {:3}\nSi/Mg_bulk = {:3}\n'.format(femg_bulk, simg_bulk)
                dat_row_header = '{0:13}{1:13}{2:13}{3:13}{4:13}'.format('mass (ME)','Radius (RE)',\
                                             'density (g/cc)', 'CRF', 'CMF')
                                        
                print 'header:\n\n'
                print dat_row_header
                np.savetxt(f_name, np.transpose([mass, rad[nn,:], rho_bulk, CRF_grid, CMF_grid]), \
                            delimiter = '    ',  fmt = '%-10.4f', header = header+dat_row_header)
            nn+=1 #count which composition  
               
                           
    
    #plot after making data files?
    
    if kwargs.get('plot') == True:
        import out
        out.pltmvr(mass = M, radius = R, labels = label)
    return
    



########################################################################
'''
grids for a range of masses and composition. Compositions are given for mantle
and core mass fraction is given a an input
Each composition is one grid
'''
########################################################################

def mvr_grid_cmf(**kwargs):
    
    
    
    #FeMg, SiMg, CMF arrays in form (start,stop, interval)
    mass = np.arange(*kwargs.get('mass'))
    FeMg = np.arange(*kwargs.get('femg'))
    SiMg = np.arange(*kwargs.get('simg'))
    CMF  = np.arange(*kwargs.get('cmf'))
    
    if kwargs.get('h2o') == None:
        wt_h2o = [0.0]
    else:
        wt_h2o = np.arange(*kwargs.get('h2o'))
    
    #filenames come in as lists or arrays

    if kwargs.get('filename') == None:
        default_filename = True
    else:
        default_filename = False
    
   
    
    CaMg       = 0.06
    AlMg       = 0.08
    Fe_ox      = 0.
    wt_Sic, wt_Oc , wt_Sc = 0.,0.,0.
    

    ##arrays to store code output for plotting and file output
    ##sample exoplex for all these cases mass_grid     = np.zeros(grid_size)
    
    #number of rows per file
    n_rows = int(len(mass))
    
    #number of files (each composition gets one)
    n_files = int(len(FeMg)*len(SiMg)*len(CMF)*len(wt_h2o))
    
    rad            = np.zeros((n_files,n_rows))
    rho_bulk       = np.zeros(n_rows)
    CRF_grid       = np.zeros(n_rows)
    CMF_grid       = np.zeros(n_rows)
    
    R     = np.empty(n_files, dtype = object)
    RR     = np.empty(n_files, dtype = object)
    M     = np.empty(n_files, dtype = object)
    R_test = np.empty(0, dtype = object)
    label = np.empty(n_files, dtype = 'S30')
    nn = 0
    
    for i in range(len(CMF)):
        for j in range(len(FeMg)):
            for k in range(len(SiMg)):
                for h in range(len(wt_h2o)):
                    
                    #set file names for this composition
                    if default_filename:
                        f_temp = 'SiMg_FeMg_CMF_h2o_{:3}_{:3}_{:3}_{:3}'.format(SiMg[k], FeMg[j],CMF[i], wt_h2o[h])
                        f_name = 'Solutions/Grids/'+f_temp.replace('.', ',')+'.dat'
                    else:
                        f_name = kwargs.get('filename') 
                  
                    #initialize exoplex inputs for this composition
                    Fix_core = {'fix_man': True, 'wtCore': CMF[i]}

                    compositional_params = [wt_h2o[h],FeMg[j],SiMg[k],CaMg,AlMg, Fe_ox ,wt_Sic, \
                                  wt_Oc , wt_Sc]

                    structure_params =  [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                                 Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                                 Core_rad_frac_guess,Mantle_potential_temp, h20_radfrac_guess, T_surface_h2o]

                    layers = [num_mantle_layers,num_core_layers,number_h2o_layers]

                    sol_filename = 'NXNN00'

                    #run exoplex for range of masses where composition is fixed above
                    for m in range(n_rows):

                        #exoplex
                        Planet = exo.run_planet_mass(mass[m], compositional_params,structure_params,layers,sol_filename, Fix_core)
                        
                        rad[nn,m]      = (Planet['radius'][-1])/REarth #R_Earth units
                        rho_bulk[m] = Planet['mass'][-1]/(4./3 *np.pi*(np.power(Planet['radius'][-1],3))) #g/cc
                        CRF_grid[m] = Planet['radius'][num_core_layers]/Planet['radius'][-1]
                        CMF_grid[m] = CMF[i]
                    
                    #store mass and radius arrays in an object array for plotting
                    M[nn]  = mass
                    R[nn]  = rad[nn,:]
                    #RR[nn] = rad
                    label[nn] = 'SiMg={:1},FeMg={:1},CMF={:1}'.format(SiMg[k], FeMg[j], CMF[i])
                    
                    #pdb.set_trace()
                    #bulk composition for file header
                    femg_bulk     = Planet['bulk_ratios'][0]
                    simg_bulk     = Planet['bulk_ratios'][1]
                        
                    header = 'Fe/Mg_bulk = {:3}\nSi/Mg_bulk = {:3}\n'.format(femg_bulk, simg_bulk)
                    dat_row_header = '{0:13}{1:13}{2:13}{3:13}{4:13}'.format('mass (ME)','Radius (RE)',\
                                                 'density (g/cc)', 'CRF', 'CMF')
                                            
                    print 'header:\n\n'
                    print dat_row_header
                    np.savetxt(f_name, np.transpose([mass, rad[nn,:], rho_bulk, CRF_grid, CMF_grid]), \
                                delimiter = '    ',  fmt = '%-10.4f', header = header+dat_row_header)
                nn+=1 #count which composition  
                   
                               
    
    #plot after making data files?
    
    if kwargs.get('plot') == True:
        import out
        out.pltmvr(mass = M, radius = R, labels = label)
    return
    
    
    #######
    
    #Finish up the script where you plot MvR after making a grid
    
    
    
    
    
    
    
    
    
    
