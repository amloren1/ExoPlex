
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.

#####
'''
part 1 :
The goal is to fix the composition of the Earth's mantle. Normally, the code
takes in bulk composition and convolves them. In this scenario, we must 
FIX MANTLE COMPOSITION AND CORE SIZE. 
- the problem with just seperating the mantle composition is that the core
size depends on bulk composition in this model
- I could just use the same bulk raitios from the McD 03 paper but,
I think we could get the mantle and core down better.

TODO:
- 

'''
#####

import os
import sys
import numpy as np
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))

import matplotlib.pyplot as plt
import ExoPlex as exo

prem_dat = '../../PREM.csv'

REarth = 6371 #kilometers

verbose = True
def prem():
    #Column 1 R,	Earth	radius	(km)
    #Column	2 Z,	Depth	(km)
    #Column	3 rho,	Density	(gm/cc)
    #Column	4 VPV,	P-wave	velocity	(km/sec)	{in	the	vertical	direction	where	anisotropic)
    #Column	5 VPH,	P-wave	velocity	(km/sec)	{in	the	horizontal	direction	where	anisotropic)
    #Column	6 VSV,	S-wave	velocity	(km/sec)	{in	the	vertical	direction	where	anisotropic)
    #Column	7 VSH,	S-wave	velocity	(km/sec)	{in	the	horizontal	direction	where	anisotropic)
    #Column	8 epsilon,	an	anisotropy	parameter	(ignore)
    #Column	9 Qmu	Quality	factor	relating	to	rigidity,	 (ignore)
    #Column	10 Qkappa.	Quality	factor	relating	to	incompressibility,	(ignore)

    dat = np.genfromtxt(prem_dat, delimiter = ',', usecols = (0,1,2,3,4,5,6))
    
    #radius where r[0] is the center
    rad         = dat[:,0]
    depth       = dat[:,1]
    rho_depth   = dat[:,2]
    rho_rad     = rho_depth[::-1]

    print rho_depth[-1]
    print depth[-1]
    return(rad, depth, rho_depth, rho_rad)
    

#function used to exchange the input molar ratio of composition to
#mass ratios for creating a Perplex file based on input MANTLE comp
def mol_to_mass_ratios(simg, femg, camg, almg, cor_wt):
    #constants, atomic masses
    mFe      = 55.845
    mMg      = 24.306
    mSi      = 28.0867
    mO       = 15.9994
    mS       = 32.0650
    mCa      = 40.078
    mAl      = 26.981
    
    A = np.array([[0, 0, simg, -1, 0, 0],
                [-1, 0, femg, 0, 0, 0],
                    [0, 0, camg, 0, -1, 0],
                    [0, 0, almg, 0, 0, -1],
                    [mFe ,mO ,mMg ,mSi ,mCa ,mAl],
                    [ 1, -1, 1, 2, 1, 1.5]])
    
    b = np.array([0, 0, 0, 0, 100., 0])

    #mol = [nFe, nO, nMg, nSi, nCa, nAl]
    mol  = np.linalg.solve(A,b)
    
    if verbose:
        print '\nCompare composition inputs with calculated outputs:'
        print 'Fe/Mg_in = {} = {} = FeMg_calc'.format(femg, mol[0]/mol[2])
        print 'Si/Mg_in = {} = {} = SiMg_calc'.format(simg, mol[3]/mol[2])
        print 'Ca/Mg_in = {} = {} = CaMg_calc'.format(camg, mol[4]/mol[2])
        print 'Al/Mg_in = {} = {} = alMg_calc\n'.format(almg, mol[5]/mol[2])
    
    M_man = 100.
    
    #these are the inputs for perplex
    feo   = round(100*(mFe+mO)*mol[0]/M_man, 6)
    mgo   = round(100*(mMg+mO)*mol[2]/M_man, 6)
    sio2  = round(100*(mSi+2*mO)*mol[3]/M_man, 6)
    cao   = round(100*(mCa+mO)*mol[4]/M_man, 6)
    al2o3 = round(100*(mAl+1.5*mO)*mol[5]/M_man, 6)
    
    wtTot = feo+mgo+sio2+cao+al2o3
    
    if verbose:
        print '\nMantle composition input for perplex:'
        print 'FeO = {}\nMgO = {}\nSiO2 = {} \nCaO = {} \nAl2O3 = {}'.format(feo, mgo, sio2, \
            cao, al2o3)
        print '\nwtTot = {}'.format(wtTot)
    
    comp_truncate = {'FeO': feo, 'SiO2': sio2, 'MgO': mgo, \
                     'CaO': cao,'Al2O3': al2o3, 'cor_wt': cor_wt}

    return(comp_truncate)



def Earth_fix_Mcor(comosition, corComp, man_only):
    
    SiMg       = comosition.get('SiMg')
    FeMg       = comosition.get('FeMg')
    CaMg       = comosition.get('CaMg')
    AlMg       = comosition.get('AlMg')
    fFeO_m     = comosition.get('fFeO') 
    
    Mass = 1.0 #Earth!


    #name of solution? might just do away with this, store star names in some 
    #directory within the solution files
    Star = 'Sun'
    
    
    
    #describe mantle compoisition independently?
    #if yes, run  routine to find inputs for perplex
    #otherwise, keep going, all will be convolved in functions.get_percents
    if man_only.get('fix_man'):
        M_frac_cor = man_only.get('wtCore')
        comp_model = mol_to_mass_ratios(SiMg, FeMg, CaMg, AlMg, M_frac_cor)
    else:
        comp_model = False

    wt_frac_Si_core = corComp.get('wtSi')
    wt_frac_O_core  = corComp.get('wtO')
    wt_frac_S_core  = corComp.get('wtS')


    wt_frac_water      = 0.00 # Earth = 0.0002
    

    #(Perplex) P&T parameter space definitions for perplex
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


    #lists of compositional and structural inputs used to build planet
    compositional_params = [wt_frac_water,FeMg,SiMg,CaMg,AlMg,fFeO_m ,wt_frac_Si_core, \
                          wt_frac_O_core,wt_frac_S_core]


    structure_params =  [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                         Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                         Core_rad_frac_guess,Mantle_potential_temp, h20_radfrac_guess, T_surface_h2o]

    layers = [num_mantle_layers,num_core_layers,number_h2o_layers]
    
    filename = Star
    
    
    #run exoplex model here
    Planet = exo.run_planet_mass(Mass, compositional_params,structure_params,layers,filename, comp_model)
    
    
    
    #setup plots
    fig, ax =  plt.subplots(figsize = (15,10))
    
    #import prem data (rad, depth, rho_depth, rho_rad)
    prem_dat = prem()
    
    depth_prm = prem_dat[1]
    rho_dep   = prem_dat[2]
        
    ax.plot(depth_prm, rho_dep, label = 'PREM',  lw = 5, ls = '-.', color = 'black')
    
    plt.rc('font', family='serif')

    
    #ax.set_ylim(0., (max(Planet['density']) / 1.e3) + 1.)
    lab_size = 23
    tic_size = 18
    ax.set_xlim(0., max(depth_prm))
    ax.set_ylabel("Density (g/cm$^3$)", fontsize = lab_size )
    ax.set_xlabel("Depth (km)", fontsize = lab_size)
    ax.tick_params(direction='in', length=6, labelsize = tic_size)
    ax.grid(color='grey', linestyle='-', alpha = 0.4, linewidth=.7)

    #data to plot
    depth_ep = Planet['radius'][-1]-Planet['radius']
    depth_ep = depth_ep/1e3
    rho_ep = Planet['density']/1e3
    data = np.array([depth_ep, Planet['density']/1e3])
    data = data.T

    #send model data to a file
    np.savetxt('Earth_McD_bulk.dat', data, fmt = '%.4f', delimiter = '\t', header = "depth (km)\tdensity (kg/m^3)")

    
    #Plot the Exoplex model
    ax.plot(depth_ep, rho_ep, label = 'ExoPlex', lw = 4, color = 'magenta')
    
    
    print 'radius of planet'
    print Planet['radius'][-1]/1000
    
    #print this stuff to make sure you are not going insane in da membrane
    print
    print "Mass = ", '%.3f'%(Planet['mass'][-1]/5.97e24), "Earth masses"
    print "Core Mass Fraction = ", '%.3f'%(100.*Planet['mass'][num_core_layers]/Planet['mass'][-1])
    print "Core Radius Fraction = ", '%.3f'%(100.*Planet['radius'][num_core_layers]/Planet['radius'][-1])
    print "CMB Pressure = " ,'%.3f' % (Planet['pressure'][num_core_layers]/10000), "GPa"
    print "Central pressure = {} GPa".format(Planet['pressure'][0]/10000)
    path_to_figs = '/home/alejandro/Documents/M-R Stuff/ThesisFigs'

    plt.legend(loc = 'lower right', fontsize = tic_size)
    plt.savefig(path_to_figs+'/Earth_McD_Bulk_XFeO.png')
    plt.show()



    



wtCore = .323 #McDonough 03
    
#bulk earth ratios from McDonough O3. 
#ran with fFeO = 0 and fFeO =  0.1324012  
bulk_Earth_comp = {'FeMg': 2.969696 , 'SiMg': 0.90303030 , \
                 'AlMg': 0.090909090 , 'CaMg': 0.06666 , 'fFeO': 0.0}

bulk_inputs = {'fix_man': False, 'wtCore': None}


bulk_Earth_comp_fFeO = {'FeMg': 2.969696 , 'SiMg': 0.90303030 , \
                 'AlMg': 0.090909090 , 'CaMg': 0.06666 , 'fFeO': 0.13240}

#earth mantle ratios from McDonough O3
#fFeO already considered in ratios
mantle_Earth_comp = {'FeMg': 0.121212121 , 'SiMg': 0.79797979797,  \
              'AlMg': 0.09090909 , 'CaMg': 0.0656565, 'fFeO': 0.0}


#Solar values from Lodders 2003
#use Al/Mg, Ca/Mg composition for Earth mantle but no FeO
solar_comp = {'FeMg': 0.813 , 'SiMg': 0.955 , \
                'AlMg': 0.090909090 , 'CaMg': 0.0656565 , 'fFeO': 0.0}
solar_input = {'fix_man': False, 'wtCore': None}
Fe_only_core = {'wtSi': 0.0, 'wtO': 0.0, 'wtS':0.0}


#core composition from McDonough 03 no O model
light_core_composition = {'wtSi': 0.06, 'wtO': 0.0, 'wtS':0.019}
man_only = {'fix_man': True, 'wtCore': 0.323}


#Earth_fix_Mcor(bulk_Earth_comp, core_composition, man_only)

#(comosition, corComp, man_only)
#####################
#1) Run Earth with max info

#Earth_fix_Mcor(mantle_Earth_comp, light_core_composition, man_only)

#2) Earth with knowledge of its bulk composition only
#use McD values with and without Fe in mantle
#left light core just because 
Earth_fix_Mcor(bulk_Earth_comp_fFeO, Fe_only_core, bulk_inputs)


#3) Earth mass with solar composition +- errors
#Earth_fix_Mcor(solar_comp, Fe_only_core, solar_input)



##################
