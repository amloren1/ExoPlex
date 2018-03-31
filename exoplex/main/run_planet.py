# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.

import os
import sys

import numpy as np

path1 = os.path.dirname(os.path.realpath(__file__))
path = path1[0:-5]
#DEBUG
#print '\n\nthis is the run_planet current path\n'
#print path
sys.path.insert(1, path)

#sys.exit()
# hack to allow scripts to be placed in subdirectories next to ExoPlex:
if not os.path.exists('main') and os.path.exists('../main'):
    sys.path.insert(1, os.path.abspath('..'))
    
import minphys
import functions
import run_perplex
#run perplex only or run full ExoPlex model?
import pdb
import params

multi_process = params.multi_process
perplex_only  = params.perplex_only
update_minphys = params.update_minphys
import multiprocessing as mp


def check_input_consistency(compositional_params, structure_params, layers):


    num_mantle_layers, num_core_layers, number_h2o_layers = layers

    wt_frac_water, FeMg, SiMg, CaMg, AlMg, mol_frac_Fe_mantle, wt_frac_Si_core, \
     wt_frac_O_core, wt_frac_S_core = compositional_params


    #Safety check for non-matching inputs, probably just remove this
    if wt_frac_water == 0. and number_h2o_layers > 0:
       print '#============================#'
       print "\n***Build error: excess in water layers for water mass fraction:\n wt_h2o = 0 wt%***"
       print "Solution: removing top {} layers allocated to water layer".format(number_h2o_layers)
       print '#============================#'
       number_h2o_layers     = 0
       #water_thickness_guess = 0
    elif wt_frac_water > 0 and number_h2o_layers == 0:
        print '#============================#'
        print "***Build error: no layers for water mass fraction:\n wt_h2 > 0 wt%***"
        print "Solution: Adding 100 layers for water envelope"
        print '#============================#'
        number_h2o_layers = 100

    new_layers = [num_mantle_layers, num_core_layers, number_h2o_layers]

    return(new_layers)

def run_planet_radius(radius_planet, compositional_params, structure_params, layers,filename, truncate_comp):
    
    print '\n*****************************************'
    print 'Running model based off of input Radius'
    print '-------------------------------------------'

    #find compositional percentages: abun. of each element, core mass frac, core composition and Perplex inputs

    #DEBUG
    #print 'Core comp'
    #print Core_wt_per.get('Fe')
    #print Core_wt_per.get('Si')
    #print Core_wt_per.get('S')
    #print Core_wt_per.get('O')

    layers = check_input_consistency(compositional_params, structure_params, layers)


    if truncate_comp.get('fix_man') == False:
        print '\nUsing molar ratios as bulk planet composiiton\n'
        print '*********************************\n'
        Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = functions.get_percents(*compositional_params)

    else:
        
        core_mass_frac = truncate_comp.get('wtCore')

        print '\nUsing molar ratios to define mantle only'
        print 'Entered core wt%% = {}\n'.format(round(core_mass_frac*100.,6))
        print '\n*****************************************'
        Mantle_wt_per, bulk_ratios = functions.get_mantle_percents(compositional_params, core_mass_frac)
        


        #set list of mantle ratios to be consistent with bulk compostiion version

        Mantle_ratios = compositional_params[1:5]


        print 'XFeO = {}'.format(bulk_ratios[4])


        wtFe = round(100.*(1- sum(compositional_params[6:])), 8)

        Core_wt_per = {'Fe': wtFe,'Si':round(100.*compositional_params[6],8) \
          ,'O':round(100.*compositional_params[7],8),'S':round(100.*compositional_params[8],8)}
        Core_mol_per = 'place holder?'


    print 'Core_wt_per: ' , Core_wt_per
    print 'Mantle_wt_per: ', Mantle_wt_per
    print 'Core_mol_per: ', Core_mol_per
    print 'core_mass_frac: ', core_mass_frac
    
    #Run perplex either in series or parallel for upper and lower mantle
    if multi_process:
        #must generate filenames with seperate function due to IO issues
        #with multiprocessing

        #(Perplex)Run fine mesh grid, Upper mantle mineralogy
        upper_man_file = functions.solfile_name(*[Mantle_wt_per,Mantle_ratios, \
        [structure_params[0],structure_params[1],structure_params[2]],filename,True])

        lower_man_file = functions.solfile_name(*[Mantle_wt_per,Mantle_ratios, \
        [structure_params[3],structure_params[4],structure_params[5]],filename,False])

        if update_minphys:
            os.remove('Solutions/'+upper_man_file+'_UM.tab')
            os.remove('Solutions/'+lower_man_file+'_LM.tab')
        #setup and run lower and upper mantle .tab files simultaneously
        p_LM = mp.Process(target = run_perplex.run_perplex, args = ([Mantle_wt_per,lower_man_file, \
                        [structure_params[3],structure_params[4],structure_params[5]],filename,False]))

        p_UM = mp.Process(target = run_perplex.run_perplex, args=([Mantle_wt_per,upper_man_file , \
                                [structure_params[0],structure_params[1],structure_params[2]],filename,True]))
        p_UM.start()
        p_LM.start()
        p_UM.join()
        p_LM.join()

    else:
        
        upper_man_file = functions.solfile_name(*[Mantle_wt_per,Mantle_ratios, \
        [structure_params[0],structure_params[1],structure_params[2]],filename,True])

        lower_man_file = functions.solfile_name(*[Mantle_wt_per,Mantle_ratios, \
        [structure_params[3],structure_params[4],structure_params[5]],filename,False])
        
        
        if update_minphys:
            os.remove('Solutions/'+upper_man_file+'_UM.tab')
            os.remove('Solutions/'+lower_man_file+'_LM.tab')
        
        LM = run_perplex.run_perplex(*[Mantle_wt_per,lower_man_file, \
                        [structure_params[3],structure_params[4],structure_params[5]],filename,False])

        UM  = run_perplex.run_perplex(*[Mantle_wt_per,upper_man_file, \
                                [structure_params[0],structure_params[1],structure_params[2]],filename,True])


        
    #only make perplex files?
    if perplex_only:
        return

    upper_man_file = 'Solutions/'+upper_man_file
    lower_man_file = 'Solutions/'+lower_man_file
    ##store upper mantle data grids: T, P, rho etc.
    grids_low, names = functions.make_mantle_grid(upper_man_file,True)
    names.append('Fe')

    #if there is a water mass fraction 0, then append h2o phases to phase list
    if layers[-1] > 0:
        names.append('liq_water')
        names.append('ice_VII')
        names.append('ice_VI')
        names.append('ice_Ih')


    #lower mantle data grids
    grids_high = functions.make_mantle_grid(lower_man_file,False)[0]


    ###Append low and high res grids. grids are rho, alpha, Cp, T,P from perplex solution
    grids = [grids_low,grids_high]

    #find mass of planet as a function of radius and composition
    Planet = functions.find_Planet_radius(radius_planet, core_mass_frac,structure_params, compositional_params, grids, Core_wt_per, layers)


    Planet['mass'] = minphys.get_mass(Planet,layers)

    Planet['phases'] = functions.get_phases(Planet, grids, layers)

    Planet['phase_names'] = names

    Planet['Vphi'], Planet['Vp'], Planet['Vs'] = functions.get_speeds(Planet, Core_wt_per, grids, layers)
    Planet['alpha'] = functions.get_alpha(Planet, Core_wt_per, grids, layers)


    return Planet













def run_planet_mass(mass_planet, compositional_params, structure_params, layers,filename, truncate_comp):
    
    print '\n*****************************************'
    print 'Running model based off of input MASS'
    print '-------------------------------------------'

    #find compositional percentages: abun. of each element, core mass frac, core composition and Perplex inputs

    params_copy = [mass_planet, compositional_params, structure_params, layers,filename, truncate_comp]
    #DEBUG
    #print 'Core comp'
    #print Core_wt_per.get('Fe')
    #print Core_wt_per.get('Si')
    #print Core_wt_per.get('S')
    #print Core_wt_per.get('O')

    layers = check_input_consistency(compositional_params, structure_params, layers)

    if truncate_comp.get('fix_man') == False:
        print '\nUsing molar ratios as bulk planet composiiton\n'
        print '\n*****************************************'
        Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac, Mantle_ratios = functions.get_percents(*compositional_params)

        bulk_ratios = compositional_params[1:6]

    else:
        core_mass_frac = truncate_comp.get('wtCore')

        #print '\n*********************************'
        print '\nUsing molar ratios to define mantle only'
        print 'Entered core wt%% = {}\n'.format(round(core_mass_frac*100.,6))
        print '\n*****************************************'

        Mantle_wt_per, bulk_ratios = functions.get_mantle_percents(compositional_params, core_mass_frac)

        #set list of mantle ratios to be consistent with bulk compostiion version

        Mantle_ratios = compositional_params[1:5]


        print 'XFeO = {}'.format(bulk_ratios[4])


        wtFe = round(100.*(1- sum(compositional_params[6:])), 8)

        Core_wt_per = {'Fe': wtFe,'Si':round(100.*compositional_params[6],8) \
          ,'O':round(100.*compositional_params[7],8),'S':round(100.*compositional_params[8],8)}
        Core_mol_per = ''


    print 'Core_wt_per: ' , Core_wt_per
    print 'Mantle_wt_per: ', Mantle_wt_per
    print 'Core_mol_per: ', Core_mol_per
    print 'core_mass_frac: ', core_mass_frac
    
    #DEBUG
    import pdb
    #pdb.set_trace()
    #if True:
    #    return

    if multi_process:
        #must generate filenames with seperate function due to IO issues
        #with multiprocessing

        #(Perplex)Run fine mesh grid, Upper mantle mineralogy
        upper_man_file = functions.solfile_name(*[Mantle_wt_per,Mantle_ratios, \
        [structure_params[0],structure_params[1],structure_params[2]],filename,True])

        lower_man_file = functions.solfile_name(*[Mantle_wt_per,Mantle_ratios, \
        [structure_params[3],structure_params[4],structure_params[5]],filename,False])

        if update_minphys:
            os.remove('Solutions/'+upper_man_file+'_UM.tab')
            os.remove('Solutions/'+lower_man_file+'_LM.tab')
        #setup and run lower and upper mantle .tab files simultaneously
        p_LM = mp.Process(target = run_perplex.run_perplex, args = ([Mantle_wt_per,lower_man_file, \
                        [structure_params[3],structure_params[4],structure_params[5]],filename,False]))

        p_UM = mp.Process(target = run_perplex.run_perplex, args=([Mantle_wt_per,upper_man_file , \
                                [structure_params[0],structure_params[1],structure_params[2]],filename,True]))


        p_UM.start()
        p_LM.start()
        p_UM.join()
        p_LM.join()

    else:
        upper_man_file = functions.solfile_name(*[Mantle_wt_per,Mantle_ratios, \
            [structure_params[0],structure_params[1],structure_params[2]],filename,True])

        lower_man_file = functions.solfile_name(*[Mantle_wt_per,Mantle_ratios, \
            [structure_params[3],structure_params[4],structure_params[5]],filename,False])

        if update_minphys:
            os.remove('Solutions/'+upper_man_file+'_UM.tab')
            os.remove('Solutions/'+lower_man_file+'_LM.tab')
            
        LM = run_perplex.run_perplex(*[Mantle_wt_per,lower_man_file, \
                        [structure_params[3],structure_params[4],structure_params[5]],filename,False])

        UM  = run_perplex.run_perplex(*[Mantle_wt_per,upper_man_file, \
                                [structure_params[0],structure_params[1],structure_params[2]],filename,True])


    #only make perplex files?
    if perplex_only:
        return

    upper_man_file = 'Solutions/'+upper_man_file
    lower_man_file = 'Solutions/'+lower_man_file
    ##store upper mantle data grids: T, P, rho etc.
    grids_low, names = functions.make_mantle_grid(upper_man_file,True)
    
    names.append('Fe')

    #if there is a water mass fraction 0, then append h2o phases to phase list
    if layers[-1] > 0:
        names.append('liq_water')
        names.append('ice_VII')
        names.append('ice_VI')
        names.append('ice_Ih')


    #lower mantle data grids
    grids_high = functions.make_mantle_grid(lower_man_file,False)[0]


    ###Append low and high res grids. grids are rho, alpha, Cp, T,P from perplex solution
    grids = [grids_low,grids_high]

    #find RADIUS of planet as a function of MASS and composition
    #Planet = functions.find_Planet_radius(radius_planet, core_mass_frac,structure_params, compositional_params, grids, Core_wt_per, layers)

    try:
    #(mass_planet, compositional_params, structure_params, layers,filename, truncate_comp)
        Planet = functions.R_of_M(mass_planet, core_mass_frac, structure_params, compositional_params, grids, Core_wt_per, layers)

    except IndexError:
        import time

        print '***************\nDeleting the UM solution file and starting over\n***************\n'
        print("3")
        time.sleep(1.5)
        print("2")
        time.sleep(2.5)
        print("1")
        time.sleep(.5)
        os.remove(upper_man_file+'_UM.tab')
        Planet = run_planet_mass(mass_planet, compositional_params, structure_params, layers,filename, truncate_comp)
        return Planet
    except ImportError:
        import time

        print '***************\nDeleting the LM solution file and starting over\n***************\n'
        print("3")
        time.sleep(1.5)
        print("2")
        time.sleep(2.5)
        print("1")
        time.sleep(.5)
        os.remove(upper_man_file+'_LM.tab')
        Planet = run_planet_mass(mass_planet, compositional_params, structure_params, layers,filename, truncate_comp)
        return Planet

    Planet['phases'] = functions.get_phases(Planet, grids, layers)

    Planet['phase_names'] = names

    Planet['Vphi'], Planet['Vp'], Planet['Vs'] = functions.get_speeds(Planet, Core_wt_per, grids, layers)

    Planet['bulk_ratios']   = bulk_ratios
    Planet['mantle_ratios'] = Mantle_ratios

    return Planet


def run_perplex_only(compositional_params, structure_params, layers, filename, truncate_comp):
    # find compositional percentages: abun. of each element, core mass frac, core composition and Perplex inputs


    # DEBUG
    # print 'Core comp'
    # print Core_wt_per.get('Fe')
    # print Core_wt_per.get('Si')
    # print Core_wt_per.get('S')
    # print Core_wt_per.get('O')

    layers = check_input_consistency(compositional_params, structure_params, layers)

    if truncate_comp == False:
        print '\n*********************************'
        print '\nUsing molar ratios as bulk planet composiiton\n'
        print '*********************************\n'
        Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = functions.get_percents(*compositional_params)

    else:
        core_mass_frac = truncate_comp.get('cor_wt')
        print '\n*********************************'
        print '\nUsing molar ratios to define mantle only'
        print 'Entered core wt%% = {}\n'.format(round(core_mass_frac * 100., 6))
        print '*********************************\n'
        wtFe = round(100. * (1 - sum(compositional_params[6:])), 8)

        Mantle_wt_per = truncate_comp
        Core_wt_per = {'Fe': wtFe, 'Si': round(100. * compositional_params[6], 8) \
            , 'O': round(100. * compositional_params[7], 8), 'S': round(100. * compositional_params[8], 8)}
        Core_mol_per = 'place holder?'

    print 'Core_wt_per: ', Core_wt_per
    print 'Mantle_wt_per: ', Mantle_wt_per
    print 'Core_mol_per: ', Core_mol_per
    print 'core_mass_frac: ', core_mass_frac

    # Run perplex either in series or parallel for upper and lower mantle
    if multi_process:
        # must generate filenames with seperate function due to IO issues
        # with multiprocessing

        # (Perplex)Run fine mesh grid, Upper mantle mineralogy
        upper_man_file = functions.solfile_name(*[Mantle_wt_per, compositional_params, \
                                                  [structure_params[0], structure_params[1], structure_params[2]],
                                                  filename, True])

        lower_man_file = functions.solfile_name(*[Mantle_wt_per, compositional_params, \
                                                  [structure_params[3], structure_params[4], structure_params[5]],
                                                  filename, False])

        # setup and run lower and upper mantle .tab files simultaneously
        p_LM = mp.Process(target=run_perplex.run_perplex, args=([Mantle_wt_per, compositional_params, \
                                                                 [structure_params[3], structure_params[4],
                                                                  structure_params[5]], filename, False]))

        p_UM = mp.Process(target=run_perplex.run_perplex, args=([Mantle_wt_per, compositional_params, \
                                                                 [structure_params[0], structure_params[1],
                                                                  structure_params[2]], filename, True]))

        p_UM.start()
        p_LM.start()
        p_UM.join()
        p_LM.join()

    else:
        lower_man_file = run_perplex.run_perplex([Mantle_wt_per, compositional_params, \
                                                  [structure_params[3], structure_params[4], structure_params[5]],
                                                  filename, False])

        upper_man_file = run_perplex.run_perplex([Mantle_wt_per, compositional_params, \
                                                  [structure_params[0], structure_params[1], structure_params[2]],
                                                  filename, True])














