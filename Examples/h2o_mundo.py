import os
import sys
import numpy as np
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))


import ExoPlex as exo


MEarth = 5.972e24 #kg
REarth = 6371   #KM
def rho_bulk(M, R):

    rho = M/((4./3)*np.pi*np.pow(R,3))

    return rho


'''
if __name__ == "__main__":
    #First user must input the ratios
    Radius_planet = 1.3


    #Star = 'Sun'

    #CaMg =0.0616595

    #SiMg =0.954993
    #AlMg = 0.0851138
    #FeMg = 0.812831

    Star = 'Rando'
    CaMg = 0.0616595
    SiMg = np.arange(0.7,1.3, 0.1)
    AlMg = 0.0851138
    FeMg = np.arange(0.7,1.3, 0.1)
    # among light elements, vary only the silicon value
    wt_frac_Si_core = np.arange(0,0.3, 0.05) #0 to 30wt% Si in core
    wt_frac_O_core = 0.
    wt_frac_S_core = 0.

    #h2o content 0 to 1.0
    wt_frac_water = 0.0

    #how much Fe in the mantle by mole. This is Fe oxidation state
    mol_frac_Fe_mantle = np.arange(0,0.3, 0.05)

    #(Perplex) P&T parameter space definitions for perplex
    #UM-upper mantle, LM-lower mantle
    Pressure_range_mantle_UM = '3000 1400000'
    Temperature_range_mantle_UM = '1400 3000'
    resolution_UM = '60 60'

    Pressure_range_mantle_LM = '1250000 6500000'
    Temperature_range_mantle_LM = '2500 5000'
    resolution_LM = '50 50'


    #layers, like concentric shells set here in each region: core, mantle, h20 envelope
    num_mantle_layers = 2000
    num_core_layers = 1000
    number_h2o_layers = 0

    #temperature at surface if no water layer. Essentially temperature below the crust
    Mantle_potential_temp = 1700.

    for i in range(len(FeMg)):
        for j in range(len(SiMg)):
            for p in range(len(wt_frac_Si_core)):
                for q in range(len(mol_frac_Fe_mantle)):

                    compositional_params = [wt_frac_water,FeMg[i],SiMg[j],CaMg,AlMg,mol_frac_Fe_mantle[q],wt_frac_Si_core[p], \
                                    wt_frac_O_core,wt_frac_S_core]

                    Core_rad_frac_guess = .54
                    structure_params =  [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                                    Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                                    Core_rad_frac_guess,Mantle_potential_temp]

                    layers = [num_mantle_layers,num_core_layers,number_h2o_layers]
                    filename = Star
                    Planet = exo.run_planet_radius(Radius_planet,compositional_params,structure_params,layers,filename)



'''
#--------------------------------------------------------------------------#
# write data files for input ranges of various scenarios
#--------------------------------------------------------------------------#


#==========================================================================#
# Fix Fe,Si/Mg and XFeO
# Vary the water content and the Si_corewt%
# XFeO can be varied but should be left at 0 or XFeO_max
#==========================================================================#


def fixed_Si_Fe_XFeO(Radius, FeMg, SiMg, XFeO):

    filename = 'R_FeMg_SiMg_XFeO_%.2f_%.2f_%.2f_%.2f.dat' % \
        (Radius,FeMg, SiMg, XFeO)

    print filename
    data_file = open(filename, 'w')
    dat_row_header = '{0:10s}{1:10s}{2:10s}{3:10s}{4:10s}{5:10s}{6:10s}{7:10s}{8:10s}{9:10s}'.format('Fe/Mg', \
                    'Si/Mg' ,'XFeO', 'Si_cwt%' , 'H2Owt%', \
                 'MEarth', 'density', 'core_wt%', 'core_rad%', 'h20_rad%')
    #'Fe/Mg    Si/Mg    XFeO    Si_cwt%    H2Owt%    MEarth    density    core_wt%    core_rad%    h20_rad%'
    data_file.write(dat_row_header)
    Radius_planet = Radius

    Star = 'Rando'
    CaMg = 0.0616595
    AlMg = 0.0851138

    # among light elements, vary only the silicon value
    wt_frac_Si_core = np.arange(0.1,0.3, 0.05) #0 to 30wt% Si in core
    wt_frac_O_core = 0.
    wt_frac_S_core = 0.

    #h2o content 0 to 1.0
    wt_frac_water = np.arange(0,0.5, 0.05)

    #how much Fe in the mantle by mole. This is Fe oxidation state
    mol_frac_Fe_mantle = XFeO

    #(Perplex) P&T parameter space definitions for perplex
    #UM-upper mantle, LM-lower mantle
    Pressure_range_mantle_UM = '3000 1400000'
    Temperature_range_mantle_UM = '1400 3000'
    resolution_UM = '60 60'

    Pressure_range_mantle_LM = '1250000 6500000'
    Temperature_range_mantle_LM = '2500 5000'
    resolution_LM = '50 50'


    #layers, like concentric shells set here in each region: core, mantle, h20 envelope
    num_mantle_layers = 2000
    num_core_layers = 1000
    number_h2o_layers = 200
    n_tot = num_mantle_layers + num_core_layers + number_h2o_layers

    #temperature at surface if no water layer. Essentially temperature below the crust
    Mantle_potential_temp = 1700.


    for q in range(len(wt_frac_water)):
        for p in range(len(wt_frac_Si_core)):

            compositional_params = [wt_frac_water[q],FeMg,SiMg,CaMg,AlMg,mol_frac_Fe_mantle,wt_frac_Si_core[p], \
                            wt_frac_O_core,wt_frac_S_core]

            Core_rad_frac_guess = .54
            structure_params =  [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                            Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                            Core_rad_frac_guess,Mantle_potential_temp]

            layers = [num_mantle_layers,num_core_layers,number_h2o_layers]
            filename = Star
            Planet = exo.run_planet_radius(Radius_planet,compositional_params,structure_params,layers,filename)

            #things to output to data file
            mas = (Planet['mass'][-1])
            rho = rho_bulk(mas,Radius*REarth)
            corwt = 100.*Planet['mass'][num_core_layers]/Planet['mass'][-1]
            cor_rad = 100.*Planet['radius'][num_core_layers]/Planet['radius'][-1]
            h2o_rad = Planet['radius'][n_tot-number_h2o_layers]/Planet['radius'][n_tot]

            #here's where we write to a file
            #'Fe/Mg    Si/Mg    XFeO    Si_cwt%    H2Owt%    MEarth    density    core_wt%    core_rad%    h20_rad%'
            data_file.write('{0:10f}{1:10f}{2:10f}{3:10f}{4:10f}{5:10s}{6:10s}{7:10s}{8:10s}{9:10s}'.format(FeMg, \
                    Si/Mg ,XFeO, wt_frac_Si_core[p] , wt_frac_water[q], \
                 mas/MEarth, rho, corwt, cor_rad, h2o_rad))
            sys.exit()


fixed_Si_Fe_XFeO(1.0,1.0,1.0,0.4)



#==========================================================================#
#==========================================================================#
#==========================================================================#
#==========================================================================#
#==========================================================================#
#==========================================================================#


















