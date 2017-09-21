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


if __name__ == "__main__":

    Star = 'Rando'
    CaMg = 0.0616595
    AlMg = 0.0851138

    # among light elements, vary only the silicon value
    wt_frac_Si_core = np.arange(0.0,0.3, 0.05) #0 to 30wt% Si in core
    wt_frac_O_core  = 0.
    wt_frac_S_core  = 0.

    #h2o content 0 to 1.0
    wt_frac_water = 0.0

    #how much Fe in the mantle by mole. This is Fe oxidation state
    mol_frac_Fe_mantle = [0.0, 0.05, 0.1, 0.15, 0.20]

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

'''

    #lists of compositional and structural inputs used to build planet

    structure_params =  [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                         Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                         Core_rad_frac_guess,Mantle_potential_temp, h20_radfrac_guess, T_surface_h2o]


    for i in range(len(FeMg)):
        for j in range(len(SiMg)):
            for p in range(len(wt_frac_Si_core)):
                for q in range(len(mol_frac_Fe_mantle)):

                    compositional_params = [wt_frac_water,FeMg[i],SiMg[j],CaMg,AlMg,mol_frac_Fe_mantle[q],wt_frac_Si_core[p], \
                                    wt_frac_O_core,wt_frac_S_core]


                    layers = [num_mantle_layers,num_core_layers,number_h2o_layers]
                    filename = Star
                    Planet = exo.run_planet_radius(Radius_planet,compositional_params,structure_params,layers,filename)

'''



def rho_bulk(M, R):

    rho = M/((4./3)*np.pi*np.power(R,3))

    return round(rho,6)

#--------------------------------------------------------------------------#
# Contour plot function
#--------------------------------------------------------------------------#

def cohntor(xlab , ylab, X, Y, Z1, Z2):

    fig, (ax1, ax2) = plt.subplots(2, sharex=True ,figsize=(30,50))
    cp1 = ax1.contourf(X, Y, Z1, size = 18)
    cp2 = ax2.contourf(X, Y, Z2, size = 22)
    fig.colorbar(cp1, ax = ax1).set_label(r'Mass (M$_\oplus$)', size = 18)
    fig.colorbar(cp2, ax = ax2).set_label(r'$\rho$ (kg m$^{-3}$)', size = 18)

    ax2.set_ylim(0,max(Y[:,0]))
    ax1.set_ylim(0,max(Y[:,0]))

    ax2.tick_params(axis = 'both', size = 16,labelsize = 20)
    ax1.tick_params(axis = 'both', labelsize = 20)
    ax2.set_xlabel(xlab, fontsize = 22)
    ax1.set_ylabel(ylab, fontsize = 22)
    ax2.set_ylabel(ylab, fontsize = 22)
    plt.show()










#--------------------------------------------------------------------------#
# write data files for input ranges of various scenarios
#--------------------------------------------------------------------------#






#==========================================================================#
# Fix Fe/Mg,Si/Mg, h2owt (no water)
# Vary XFeO and the Si_corewt%
# XFeO can be varied but should be left at 0 or XFeO_max
#==========================================================================#

def fixed_Si_Fe_Noh2o(Radius, FeMg, SiMg):

    #scale up to REarth
    Radius_planet = Radius

    #labels for contour plot
    Xlab = r'Si$_{core}$ wt%'
    Ylab = r'X$_{FeO}$'

    #earth abuns from McDonough 03
    CaMg = 0.0616595
    AlMg = 0.0851138

    filename = 'R_FeMg_SiMg_CaMg_AlMg_NoH2o_{:.2f}_{:.2f}_{:.2f}_{:.2f}_{:.2f}.dat'.format(Radius,FeMg,SiMg,CaMg,AlMg)

    print filename
    data_file = open(filename, 'w')
    dat_row_header = '{0:11}{1:11}{2:11}{3:11}{4:11}{5:11}{6:11}{7:11}{8:11}{9:11}{10:11}'.format('XFeO', 'Si_cwt%','Fe/Mg', \
                    'Si/Mg', 'H2Owt%', 'MEarth', 'density', 'core_wt%', 'core_rad%', 'h20_rad%', 'calc_h20')

    data_file.write(dat_row_header)
    #vary these two

    wt_frac_Si_core = [0.0,0.01,0.02,0.05, 0.1,0.12, 0.15, 0.20, 0.25]
    XFeO            = [0.0, 0.01, 0.05,0.07,0.1,0.15,0.20,0.25]

    #make sure water parameters are zeroed out
    num_h2o_layers      = 0
    structure_params[8] = 0

    #layers in model
    layers = [num_mantle_layers,num_core_layers,num_h2o_layers]
    n_tot = sum(layers)-1

    #store independent values in arrays for contour plots
    Sic_X  = np.zeros((len(XFeO), len(wt_frac_Si_core)))
    XFeO_Y = np.zeros((len(XFeO), len(wt_frac_Si_core)))
    mas_Z  = np.zeros((len(XFeO), len(wt_frac_Si_core)))
    rho_Z  = np.zeros((len(XFeO), len(wt_frac_Si_core)))

    for q in range(len(XFeO)):
        for p in range(len(wt_frac_Si_core)):

            #store independent values in arrays for contour plots
            Sic_X[q,p]  = wt_frac_Si_core[p]
            XFeO_Y[q,p] = XFeO[q]

            compositional_params = [wt_frac_water, FeMg, SiMg, CaMg, AlMg, XFeO[q], wt_frac_Si_core[p], \
                                    wt_frac_O_core, wt_frac_S_core]

            print '\nInput composition:'
            print'\n{} = {:.2f}\n{} = {:.2f}\n{} = {:.2f}\n{} = {:.2f}\n{} = {:.4f}\n{} = {:.4f}\n'.format('XFeO ',XFeO[q],'fSic ', \
                    wt_frac_Si_core[p],'Fe/Mg', \
                    FeMg,'Si/Mg', SiMg, 'Ca/Mg', CaMg, 'Al/Mg', AlMg)

            layers = [num_mantle_layers,num_core_layers,num_h2o_layers]

            filename = 'Star'

            #run routines to build perplex solution and planet
            Planet = exo.run_planet_radius(Radius_planet,compositional_params,structure_params,layers,filename)

            continue

            #things to output to data file
            mas = (Planet['mass'][-1])
            rho = rho_bulk(mas,Radius*REarth)
            corwt = round(100.*Planet['mass'][num_core_layers]/Planet['mass'][-1],4)
            cor_rad = round(100.*Planet['radius'][num_core_layers]/Planet['radius'][-1],4)


            h2o_rad = round(abs(Planet['radius'][n_tot-num_h2o_layers]/Planet['radius'][-1]-1.00)*100.,4)
            mas_E = np.around(mas/MEarth, decimals = 6)


            h2owt_test = Planet['mass'][num_mantle_layers+num_core_layers-1]/Planet['mass'][-1]

            print '\n\nresulting h2owt = {}'.format(1.-h2owt_test)
            print 'wt frac water = {}'.format(wt_frac_water)

            #Z-axis for contour plots
            mas_Z[q,p] = mas/MEarth
            rho_Z[q,p] = rho

            #write to a file
            #'Fe/Mg    Si/Mg    XFeO    Si_cwt%    H2Owt%    MEarth    density    core_wt%    core_rad%    h20_rad%'
            data_file.write('\n{0:>10}{1:>10}{2:>10.2f}{3:>10.2f}{4:>10.2f}{5:>11.4f}{6:>11.4f}{7:>10.4f}{8:>10.4f}{9:>10.4f}{10:>10.4f}'.format(XFeO[q],\
                    wt_frac_Si_core[p]*100.,FeMg,SiMg, wt_frac_water, mas_E, rho, corwt, cor_rad, h2o_rad, (1-h2owt_test)))


    data_file.close()
    cohntor(Xlab , Ylab, Sic_X, XFeO_Y, mas_Z, rho_Z)
    return


#(Radius, FeMg, SiMg)
#fixed_Si_Fe_Noh2o(1.4, 0.5, 1.0)
fixed_Si_Fe_Noh2o(1.5, 1.5, 1.0)
sys.exit()









#==========================================================================#
# Fix Fe,Si/Mg and XFeO
# Vary the water content and the Si_corewt%
# XFeO can be varied but should be left at 0 or XFeO_max
#==========================================================================#


def fixed_Si_Fe_XFeO(Radius, FeMg, SiMg, XFeO):

    Xlab = r'Si$_{core}$ wt%'
    Ylab = r'H$_2$O wt%'

    filename = 'R_FeMg_SiMg_XFeO_%.2f_%.2f_%.2f_%.2f.dat' % \
        (Radius,FeMg, SiMg, XFeO)

    print filename
    data_file = open(filename, 'w')
    dat_row_header = '{0:11}{1:11}{2:11}{3:11}{4:11}{5:11}{6:11}{7:11}{8:11}{9:11}{10:11}'.format('Fe/Mg', \
                    'Si/Mg' ,'XFeO', 'Si_cwt%' , 'H2Owt%', \
                 'MEarth', 'density', 'core_wt%', 'core_rad%', 'h20_rad%', 'calc_h20')
    #'Fe/Mg    Si/Mg    XFeO    Si_cwt%    H2Owt%    MEarth    density    core_wt%    core_rad%    h20_rad%'
    data_file.write(dat_row_header)
    Radius_planet = Radius

    Star = 'Rando'
    #CaMg = 0.0616595
    #AlMg = 0.0851138

    CaMg = 0
    AlMg = 0

    # among light elements, vary only the silicon value
    wt_frac_Si_core = [0.0, 0.05, 0.1, 0.15, 0.20, 0.25] #0 to 30wt% Si in core
    wt_frac_Si_core = [0.15]
    wt_frac_O_core = 0.
    wt_frac_S_core = 0.

    #h2o content 0 to 1.0
    wt_frac_water = np.arange(0.0,0.5, 0.05)
    wt_frac_water = [0.05]
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


    #layers, like concentric shells set here in each region: core, mantle, h20 envelope (below)
    num_mantle_layers = 2000
    num_core_layers = 1000


    #temperature at surface if no water layer. Essentially temperature below the crust
    Mantle_potential_temp = 1700.

    #h2o potential Temp, surface temperature if there exists an h2o layer
    T_surface_h2o = 300. # K

    #initialize planet with these guesses for radial fraction of core and water layer
    Core_rad_frac_guess = .54
    h20_radfrac_guess = 0.1


    structure_params =  [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                         Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                         Core_rad_frac_guess,Mantle_potential_temp, h20_radfrac_guess, T_surface_h2o]

    counter = 1

    #store independent values in arrays for contour plots
    SicX = np.zeros((len(wt_frac_water), len(wt_frac_Si_core)))
    h2oY = np.zeros((len(wt_frac_water), len(wt_frac_Si_core)))
    masZ = np.zeros((len(wt_frac_water), len(wt_frac_Si_core)))
    rhoZ = np.zeros((len(wt_frac_water), len(wt_frac_Si_core)))

    for q in range(len(wt_frac_water)):
        for p in range(len(wt_frac_Si_core)):

            #make sure to add or zero out water layers and wrf guess
            if wt_frac_water[q] > 0.0:
                num_h2o_layers = 500
            else:
                num_h2o_layers = 0
                structure_params[8] = 0


            #store independent values in arrays for contour plots
            SicX[q,p] = wt_frac_Si_core[p]
            h2oY[q,p] = wt_frac_water[q]


            layers = [num_mantle_layers,num_core_layers,num_h2o_layers]
            n_tot = sum(layers)-1

            #composition varies in water and core light element abundance
            compositional_params = [wt_frac_water[q],FeMg,SiMg,CaMg,AlMg,mol_frac_Fe_mantle,wt_frac_Si_core[p], \
                                    wt_frac_O_core,wt_frac_S_core]

            layers = [num_mantle_layers,num_core_layers,num_h2o_layers]

            filename = Star

            #run routines to build perplex solution and planet
            Planet = exo.run_planet_radius(Radius_planet,compositional_params,structure_params,layers,filename)


            #things to output to data file
            mas = (Planet['mass'][-1])
            rho = rho_bulk(mas,Radius*REarth)
            corwt = round(100.*Planet['mass'][num_core_layers]/Planet['mass'][-1],4)
            cor_rad = round(100.*Planet['radius'][num_core_layers]/Planet['radius'][-1],4)


            h2o_rad = round(abs(Planet['radius'][n_tot-num_h2o_layers]/Planet['radius'][-1]-1.00)*100.,4)
            mas_E = np.around(mas/MEarth, decimals = 6)


            h2owt_test = Planet['mass'][num_mantle_layers+num_core_layers-1]/Planet['mass'][-1]

            print '\n\nresulting h2owt = {}'.format(1.-h2owt_test)
            print 'wt frac water = {}'.format(wt_frac_water[q])
            #Z-axis for contour plots
            masZ[q,p] = mas/MEarth
            rhoZ[q,p] = rho

            #write to a file
            #'Fe/Mg    Si/Mg    XFeO    Si_cwt%    H2Owt%    MEarth    density    core_wt%    core_rad%    h20_rad%'
            data_file.write('\n{0:>10}{1:>10}{2:>10.2f}{3:>10.2f}{4:>10.2f}{5:>11.4f}{6:>11.4f}{7:>10.4f}{8:>10.4f}{9:>10.4f}{10:>10.4f}'.format(FeMg,\
                    SiMg ,XFeO, wt_frac_Si_core[p]*100. , wt_frac_water[q]*100. , mas_E, rho, corwt, cor_rad, h2o_rad, (1-h2owt_test)))

            #make contour plots


            ##
            #
            ##

            counter += 1
            if counter < 1:

                #data_file.close()
                #cohntor(Xlab , Ylab, SicX, h2oY, masZ, rhoZ)
                break

    #close data file and plot data
    data_file.close()
    cohntor(Xlab , Ylab, SicX, h2oY, masZ, rhoZ)

#(Radius, FeMg, SiMg, XFeO)
fixed_Si_Fe_XFeO(1.0,1.0,1.0,0.0)
#fixed_Si_Fe_XFeO(1.0,1.3,1.0,0.0)


#==========================================================================#






#==========================================================================#
#==========================================================================#
#==========================================================================#
#==========================================================================#
#==========================================================================#


















