import numpy as np
import sys
from scipy import interpolate
import minphys
import os
import pdb


def solfile_name(*args):

    Mantle_wt_per = args[0]

    FeMg = args[1][0]
    SiMg = args[1][1]
    CaMg = args[1][2]
    AlMg = args[1][3]

    filename = args[3]
    UMLM     = args[4]
    
    plxMan = str(Mantle_wt_per.get('MgO')) + ' ' + str(Mantle_wt_per.get('SiO2')) + ' ' \
             + str(Mantle_wt_per.get('FeO')) + ' ' + str(Mantle_wt_per.get('CaO')) \
             + ' ' + str(Mantle_wt_per.get('Al2O3'))+ ' ' + str(0.) #last value included for Na


    #this filename convention organizes solution files by their MANTLE composition (12.05.2017)

    solfileparamsString0 = '_' + str(round(SiMg, 4)) + '_' + str(round(FeMg, 4)) + '_' + str(
        round(CaMg, 4)) + '_' + str(round(AlMg, 4))

    #old filenames
    #solfileparamsString0 = '_' + str(round(SiMg, 3)) + '_' + str(round(FeMg, 3)) + '_' + str(
    #    round(CaMg, 3)) + '_' + str(round(AlMg, 3)) \
    #                       + '_' + str(round(mol_frac_Fe_mantle, 3)) + '_' + str(round(wt_frac_Si_core, 3))


    # changes periods to commas
    solfileparamsString = solfileparamsString0.replace('.', ',')
    solutionFileNameMan = 'SiMg_FeMg_CaMg_AlMg' + solfileparamsString + '_MANTLE'

    filename = solutionFileNameMan

    #pdb.set_trace()
    return filename




verbose = True


#function used to exchange the input molar ratio of composition to
#mass ratios for creating a Perplex file based on input MANTLE comp
#this will bypass internal function which
def get_mantle_percents(args, cor_wt):
    #constants, atomic masses
    mFe      = 55.845
    mMg      = 24.306
    mSi      = 28.0867
    mO       = 15.9994
    mS       = 32.0650
    mCa      = 40.078
    mAl      = 26.981

    femg = args[1]
    simg = args[2]
    camg = args[3]
    almg = args[4]

    WtSi_c = args[6]
    WtO_c  = args[7]
    WtS_c  = args[8]

    M_tot = 100.
    M_man = (1-cor_wt)*M_tot
    M_core = cor_wt*M_tot

    A = np.array([[0, 0, simg, -1, 0, 0],
                [-1, 0, femg, 0, 0, 0],
                    [0, 0, camg, 0, -1, 0],
                    [0, 0, almg, 0, 0, -1],
                    [mFe ,mO ,mMg ,mSi ,mCa ,mAl],
                    [ 1, -1, 1, 2, 1, 1.5]])

    b = np.array([0, 0, 0, 0, M_man, 0])

    #mol = [nFe, nO, nMg, nSi, nCa, nAl]
    mol  = np.linalg.solve(A,b)


    if verbose:
        print '\nCompare composition inputs with calculated outputs:'
        print 'Fe/Mg_in = {} = {} = FeMg_calc'.format(femg, mol[0]/mol[2])
        print 'Si/Mg_in = {} = {} = SiMg_calc'.format(simg, mol[3]/mol[2])
        print 'Ca/Mg_in = {} = {} = CaMg_calc'.format(camg, mol[4]/mol[2])
        print 'Al/Mg_in = {} = {} = alMg_calc\n'.format(almg, mol[5]/mol[2])




    #these are the inputs for perplex
    feo   = round(100*(mFe+mO)*mol[0]/M_man, 6)
    mgo   = round(100*(mMg+mO)*mol[2]/M_man, 6)
    sio2  = round(100*(mSi+2*mO)*mol[3]/M_man, 6)
    cao   = round(100*(mCa+mO)*mol[4]/M_man, 6)
    al2o3 = round(100*(mAl+1.5*mO)*mol[5]/M_man, 6)

    wtTot = feo+mgo+sio2+cao+al2o3


    #calculate XFeO
    M_Fe_core      = (1.- WtSi_c- WtS_c- WtO_c)*M_core
    M_Si_core      = (WtSi_c)*M_core

    Fe_core_mol    = M_Fe_core/mFe
    Si_core_mol    = M_Si_core/mSi

    f_FeO_man = feo/100.
    Fe_man_mol  = (M_man*f_FeO_man)/(mFe+mO)

    XFeO = Fe_man_mol/(Fe_core_mol+Fe_man_mol)

    #update ratios to reflect bulk composition
    femg_blk = (mol[0]+Fe_core_mol)/(mol[2])
    simg_blk = (mol[3]+Si_core_mol)/(mol[2])


    if verbose:
        print '\nMantle composition input for perplex:'
        print 'FeO = {}\nMgO = {}\nSiO2 = {} \nCaO = {} \nAl2O3 = {}'.format(feo, mgo, sio2, \
            cao, al2o3)
        print '\nwtTot = {}'.format(wtTot)
        print '\nCalculated bulk elemental ratios'
        print 'Si/Mg = {} \nFe/Mg = {}\n'.format(simg_blk, femg_blk)



    comp_truncate = {'FeO': feo, 'SiO2': sio2, 'MgO': mgo, \
                     'CaO': cao,'Al2O3': al2o3, 'cor_wt': cor_wt}

    #update_comp = {'SiMg' : simg_blk, 'FeMg' : femg_blk, 'XFeO' : XFeO}

    bulk_ratios = [femg_blk, simg_blk, camg, almg, XFeO]
    #DEBUG
    import pdb
    #pdb.set_trace()

    #sys.exit()

    return(comp_truncate, bulk_ratios)





def get_percents(*args):
    FeMg = args[1]
    SiMg = args[2]
    CaMg = args[3]
    AlMg = args[4]
    mol_frac_Fe_mantle = args[5]
    wt_frac_Si_core    = args[6]
    wt_frac_O_core     = args[7]
    wt_frac_S_core     = args[8]

    MgSi = 1./SiMg
    FeSi = FeMg*MgSi
    CaSi = CaMg*MgSi
    AlSi = AlMg*MgSi

#constants, atomic masses
    mFe      = 55.845
    mMg      = 24.306
    mSi      = 28.0867
    mO       = 15.9994
    mS       = 32.0650
    mCa      = 40.078
    mAl      = 26.981

    # system of equation which when solved gives molar values of elements in mantle and core
    # x = [nFec,nSic,nOc, nSc | ,nFem,nMgm,nSim,nOm, nCam, nAlm]

    #Sum of all masses = 100 g
    b = np.array([0., 0. , 0. , 0. ,  0. , 0. , 0 , 0.,0., 100.])

############################################################################################

    #Comment here
    #Double checks
    #Return only moles of elements, no molecules
    A = np.array([ [0., 0., 0. , 0. , -1. ,-1., -2. , 1. , -1., -1.5] ,
            [1. , -1.*FeSi , 0. , 0. , 1. , 0. , -1.*FeSi , 0.,0.,0.] ,
           [0., -1.*MgSi, 0., 0. ,0. ,1. ,-1.*MgSi ,0.,0.,0.] ,
            [0., -1.*CaSi, 0., 0. ,0. ,0. ,-1.*CaSi ,0., 1. ,0.] ,
        [0., -1.*AlSi, 0., 0. ,0. ,0. ,-1.*AlSi ,0.,0., 1.] ,
          [mol_frac_Fe_mantle , 0. , 0. ,0., (mol_frac_Fe_mantle-1.) ,  0. , 0. , 0. ,0., 0.] ,
          [wt_frac_Si_core*mFe , (wt_frac_Si_core-1.)*mSi , wt_frac_Si_core*mO , wt_frac_Si_core*mS \
                       , 0. , 0. , 0. , 0. , 0., 0.] ,
          [wt_frac_O_core*mFe , (wt_frac_O_core)*mSi , (wt_frac_O_core-1.)*mO , wt_frac_O_core*mS
                       , 0. , 0. , 0. , 0. , 0., 0.]  ,
          [wt_frac_S_core*mFe , (wt_frac_S_core)*mSi , (wt_frac_S_core)*mO , (wt_frac_S_core-1.)*mS
                       , 0. , 0. , 0. , 0. , 0., 0.] ,
          [mFe , mSi , mO, mS ,  mFe , mMg , mSi ,mO, mCa, mAl]])


    #returns number of moles of each element between core and mantle.

    #returns number of moles assuming 100 g planet
    #

    Num_moles = np.linalg.solve(A,b)


    #To do: Adjust so oxygen is last in list
    ## find masses and wt% below for perplex ##

    #Splitting up into lists
    #in order Fe, Si, O, S in core
    Core_moles = Num_moles[:4]

    #in order Fe, Mg, Si, O, Ca, Al in mantle
    Mantle_moles = Num_moles[4:]

    tot_moles_core = sum(Core_moles)

    mass_of_Core = (mFe*(Core_moles[0])+(mSi*Core_moles[1])\
                    +((mO)*Core_moles[2])+(mS*Core_moles[3]))
    mass_of_Mantle = (mFe*Mantle_moles[0])+(mMg*Mantle_moles[1])\
                   +(mSi*Mantle_moles[2])+(mO*Mantle_moles[3])\
                   +(mCa*Mantle_moles[4])+(mAl*Mantle_moles[5])

    #print calculated molar values to screen

    if verbose:
        verbosity(Num_moles, mass_of_Core, mSi)



    Mtot= mass_of_Core+mass_of_Mantle #in g

    core_mass_frac = mass_of_Core/Mtot
    #Weight percents of mantle oxides
    #Weight percents assuming FeO, MgO, SiO2, CaO, Al2O3

    FeO_mant_wt = Mantle_moles[0]*(mFe+mO)/mass_of_Mantle
    MgO_mant_wt = Mantle_moles[1]*(mMg+mO)/mass_of_Mantle
    SiO2_mant_wt = Mantle_moles[2]*(mSi+(2.*mO))/mass_of_Mantle
    CaO_mant_wt = Mantle_moles[4]*(mCa+mO)/mass_of_Mantle
    Al2O3_mant_wt = (Mantle_moles[5]/2.)*(2.*mAl+3.*mO)/mass_of_Mantle

    #Throw exception, not if statement
    #make inequality not, absolute if. Use machine precision
    total = float(FeO_mant_wt+MgO_mant_wt+SiO2_mant_wt+CaO_mant_wt+Al2O3_mant_wt)
    if total > 1.+(5.*np.finfo(float).eps) or total < 1.-(5.*np.finfo(float).eps):
        print total
        print '\n\n Mantle wt% don\'t add to 1'
        sys.exit()

    #same as above edited for printing
    #Fix for PerPlex to input and round to the 8th decimal point
    #converts to percentages

    FeO_mant_wt = abs(round(FeO_mant_wt * 100., 8))
    SiO2_mant_wt = abs(round(SiO2_mant_wt * 100., 8))
    MgO_mant_wt  = abs(round(MgO_mant_wt * 100., 8))
    CaO_mant_wt  = abs(round(CaO_mant_wt * 100., 8))
    Al2O3_mant_wt  = abs(round(Al2O3_mant_wt * 100., 8))

    Mantle_wt_per = {'FeO': FeO_mant_wt, 'SiO2': SiO2_mant_wt, 'MgO': MgO_mant_wt, \
                     'CaO': CaO_mant_wt,'Al2O3':Al2O3_mant_wt}

    #decimal fraction of materials in CORE by mass, these are perplex inputs (hence the rounding)
    # this is the wt% of Fe and the Fe in FeSi *** so total Fe in core

    Fe_core_wt   = (Core_moles[0])*(mFe)/mass_of_Core
    Si_core_wt   = Core_moles[1]*(mSi)/mass_of_Core
    O_core_wt   = Core_moles[2]*(mO)/mass_of_Core
    S_core_wt   = Core_moles[3]*(mS)/mass_of_Core

    #Throw exception, not if statement
    #make inequality not, absolute if. Use machine precision
    corwt_tot = round(S_core_wt+O_core_wt+Si_core_wt+Fe_core_wt,8)
    if corwt_tot != 1. and corwt_tot > 1:
        print '\n\n*****Exiting program*****'
        print 'Core wt%% don\'t add up'
        print 'S_core_wt + O_core_wt + Si_core_wt + Fe_core_wt = %r' % (corwt_tot)
        print 'Siwt = %.5f' % (Si_core_wt)
        print 'Siwt_input = %.5f'% wt_frac_Si_core
        print '*************************'
        sys.exit()

    Fe_core_wt = abs(round(Fe_core_wt*100.,8))
    Si_core_wt  = abs(round(Si_core_wt*100.,8))
    O_core_wt = abs(round(O_core_wt*100.,8))
    S_core_wt  = abs(round(S_core_wt*100.,8))

    Core_wt_per = {'Fe':Fe_core_wt,'Si':Si_core_wt,'O':O_core_wt,'S':S_core_wt}
    Core_mol_per ={'Fe':Core_moles[0]/tot_moles_core,'Si':Core_moles[1]/tot_moles_core,\
                  'O':Core_moles[2]/tot_moles_core,'S':Core_moles[3]/tot_moles_core}

    #in order Fe, Mg, Si, O, Ca, Al in mantle
    Mantle_moles = Num_moles[4:]

    Mantle_ratios = [(Mantle_moles[0]/Mantle_moles[1]), (Mantle_moles[2]/Mantle_moles[1]), \
                     (Mantle_moles[4]/Mantle_moles[1]), \
                     (Mantle_moles[5]/Mantle_moles[1])]


    import pdb
    #pdb.set_trace()
    return(Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac, Mantle_ratios)







def verbosity(x, mcor, mSi):

    # MgO,SiO2,FeO, CaO, Al2O3 wt%


#    solutionFileNameMan    = 'SiMg_FeMg_CaMg_AlMg_XFeO_fSic_fOc_fSc' + solfileparamsString+'_MANTLE'
#    datafileName         = 'Mass_MgSi_FeSi_XFeO_fSic_fOc_' +repr(round(PlanetMass/MEarth,3))+solfileparamsString

#    print 'Mantle solution file: \n'
#    print solutionFileNameMan

    #sys.exit()
    #for creating data files in perplex, plcCor is an entry into the build code
    #plxCor  = repr(Fecwt)+' '+repr(O2cwt)+' '+repr(Sicwt) + ' ' + repr(S2cwt)
    #solfileparamsString0 = '_'+repr(round(fSic,3)) +'_'+repr(round(fOc,3)) + '_'+ repr(round(fSc,3))
    #solfileparamsString  = solfileparamsString0.replace('.',',')
#    datafileName         = 'Mass_MgSi_FeSi_XFeO_fSic_fOc_' +repr(normMass)+solfileparamsString
    #solutionFileNameCor     = 'fSic_fOc_fSc' + solfileparamsString+'_CORE'

     #print some stuff here
    #if verbose:
    ##    print 'Mantle composition (wt%%): \n'
     #   print 'MgO = %.4f \tSiO2 = %.4f \tFeO = %.4f \tCaO = %.4f \tAl2O3 = %.4f\n' %                        \
     #           (np.abs(MgOmwt),np.abs(SiO2mwt),np.abs(FeOmwt),np.abs(CaOmwt), np.abs(Al2O3mwt) )

      #  print 'Core composition (wt%%): \n'
      #  print 'O2 = %.4f \tSi = %.4f \tS2 = %.4f \tFe = %.4f (Fe alone & Fe in FeSi)\n' %          \
       #      (np.abs(O2cwt),np.abs(Sicwt),np.abs(S2cwt),np.abs(Fecwt))

        #print 'wt%% of Core: %r \t wt%% of Mantle: %r' % (round(masfCor*100,5), round(masfMan*100,5))
        #print '*---------------------------------------------------*\n'



    # x = [nFec,nMgc,nSic,nOc, nSc | ,nFem,nMgm,nSim,nOm]
    #DEBUG:
    #DEBUG: CHECK TOTAL WEIGHT PERCENTAGE = 100%
    #print 'total mass of planet = %r' % ((masCor+masMan)/Mp)
    #print 'total wt%% in mantle = %.30f' % (FeOmwt+SiO2wt+MgOwt)
    #print 'total wt%% in core = %.30f' % ((O2cwt+Sicwt+Fecwt))

    #testValue = x[2]*(mSi+mFe)/(mFe*x[0]+mSi*x[2]+mO*x[3])
    #testValue1 = x[3]*(mO+mFe)/(mFe*x[0]+mSi*x[2]+mO*x[3])
    #print testValue1
    #print testValue

    #DEBUG: print out the ratios of Mg/Si and Fe/Si
    # x = [nFec,nSic,nOc, nSc | ,nFem,nMgm,nSim,nOm, nCam, nAlm]
    if verbose:
        print 'Calculated composition:'
        print 'XFeO  = %.4f' % (x[4]/(x[0]+x[4]))
        print 'fSic  = %.4f' % ((x[1]*(mSi))/mcor)
        print 'Si/Mg = %.4f' % ((x[1]+x[6])/x[5])
        print 'Fe/Mg = %.4f' % ((x[0]+x[4])/x[5])
        print 'Ca/Mg = %.4f' % ((x[8])/x[5])
        print 'Al/Mg = %.4f' % ((x[9])/x[5])

        #print 'fSic = %.4f' % ((x[1]*(mSi))/masfCor)
        #print 'fOc = %.4f' % ((x[2]*(mO))/masfCor)
        #print 'fSc = %.4f' % ((x[3]*(mS))/masfCor)
      #  print 'Core wt%% of Si = %.4f %%' % (100*(x[1]*(mSi))/masfCor)
        #print 'Core wt%% of O = %.4f %%' % (100*(x[2]*(mO))/masfCor)
        #print 'Core wt%% of S = %.4f %%' % (100*(x[3]*(mO))/masfCor)
        #print 'Total light element wt%% in core: %.4f %%' % (100*(x[2]*(mO)+x[1]*mSi+x[3]*mS)/masfCor)
        #print 'Mantle Mg/Si = %.4f' % (x[5]/x[6])
        #print '\nEntered values:'
        #print 'Si/Mg = %r' % SiMg
        #print 'Fe/Mg = %r' % FeMg
        #print 'Ca/Mg = %r' % CaMg
        #print 'Al/Mg = %r' % AlMg
        #print 'XFeO = %.4f' % XFeO
        #print 'fSic = %.4f' % (fSic)
        #print 'fOc = %.4f' % fOc
        #print 'fSc = %.4f' % fSc

    #sys.exit()
#  return values: masfCor,masfMan, MgOmwt,SiO2mwt,FeOmwt, Sicwt, Fecwt,O2cwt
    #must rememer: the masfMan and Cor are wt%s
    return

def make_mantle_grid(Mantle_filename,LMUM):

    if LMUM == True:
        
        file = open(Mantle_filename+'_UM.tab','r')
    else:
        file = open(Mantle_filename+'_LM.tab','r')

    temp_file = file.readlines()
    num_rows = len(temp_file[13:])
    num_columns = len(temp_file[12].split())

    header = temp_file[12].strip('\n').split()
    Phases = header[8:]

    for i in range(len(Phases)):
        Phases[i] = Phases[i].strip(",mo%")

    data = temp_file[13:]
    grid = np.zeros((num_rows,num_columns))

    for i in range(num_rows):
        #for j in range(num_columns):
        columns = data[i].strip('\n').split()
        grid[i] = [float(j) for j in columns]


    num_phases = len(grid[0][8:])
    phases_grid = np.zeros((num_rows,num_phases))
    for i in range(num_rows):
        phases_grid[i] = grid[i][8:]

    temperature_grid = [row[0] for row in grid]
    pressure_grid = [row[1] for row in grid]
    density_grid = [row[2] for row in grid]
    speed_grid = [[row[3],row[4],row[5]] for row in grid]
    alpha_grid = [row[6] for row in grid]
    cp_grid = [row[7] for row in grid]
    phase_grid = [row[8:] for row in grid]

    keys = ['temperature','pressure','density','speeds','alpha','cp','phases']
    return dict(zip(keys,[temperature_grid,pressure_grid,density_grid,speed_grid,alpha_grid,cp_grid,phase_grid])),Phases

def get_phases(Planet,grids,layers):

    num_mantle_layers, num_core_layers, number_h2o_layers = layers


    mantle_pressures = Planet['pressure'][num_core_layers:(num_core_layers+num_mantle_layers)]
    mantle_temperatures = Planet['temperature'][num_core_layers:(num_core_layers+num_mantle_layers)]
    core_pressures = Planet['pressure']
    core_temperatures = Planet['temperature']

    P_points_UM = []
    T_points_UM = []

    for i in range(len(mantle_pressures)):
        if mantle_pressures[i]<=1250000.:
            P_points_UM.append(mantle_pressures[i])
            T_points_UM.append(mantle_temperatures[i])



    Mantle_phases_UM = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']), grids[0]['phases'],
                         (P_points_UM, T_points_UM), method='linear')

    Mantle_phases_LM = [grids[0]['phases'][-1] for i in range(num_mantle_layers-len(P_points_UM))]

    if (num_mantle_layers-len(P_points_UM)) > 0:
        Mantle_phases = np.concatenate((Mantle_phases_LM,Mantle_phases_UM),axis=0)

    else:
        Mantle_phases = Mantle_phases_UM

    for i in range(len(Mantle_phases)):
        entry = np.zeros(len(Mantle_phases[i]))
        tot = sum(Mantle_phases[i])
        for j in range(len(Mantle_phases[i])):
            entry[j] = 100.*Mantle_phases[i][j]/tot
        Mantle_phases[i] = entry

    if number_h2o_layers>0:
        Phases = np.zeros((len(Planet['pressure']),len(Mantle_phases[0])+5))
    else:
        Phases = np.zeros((len(Planet['pressure']),len(Mantle_phases[0])+1))

    Core_phases = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']), grids[0]['phases'],
                         ([0 for i in range(num_core_layers)], [0 for i in range(num_core_layers)]), method='linear')




    if number_h2o_layers >0:
        Water_phases = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']), grids[0]['phases'],
                         ([0 for i in range(number_h2o_layers)], [0 for i in range(number_h2o_layers)]), method='linear')
        Phases[0:num_core_layers,0:len(Mantle_phases[0])] = Core_phases
        Phases[num_core_layers:(num_mantle_layers+num_core_layers),0:len(Mantle_phases[0])] = Mantle_phases
        Phases[(num_mantle_layers+num_core_layers):num_core_layers+num_mantle_layers+number_h2o_layers,0:len(Mantle_phases[0])] = Water_phases

    else:
        Phases[0:num_core_layers, 0:len(Mantle_phases[0])] = Core_phases
        Phases[num_core_layers:(num_mantle_layers + num_core_layers), 0:len(Mantle_phases[0])] = Mantle_phases

    # add in the core:
    if number_h2o_layers > 0:
        for i in range(len(Phases)):
            if i < num_core_layers-1:
                #append Fe, liq_water, ice_vii, ice_vi, ice_ih
                #np.resize(Phases[i],len(Phases[i])+4)
                Phases[i][len(Mantle_phases[0]):] =[100, float('NaN'), float('NaN'), float('NaN'),float('NaN')]
            elif i < (num_mantle_layers+num_core_layers):
                 #Phases[i] =np.resize(Phases[i],len(Phases[i])+4)

                 Phases[i][len(Mantle_phases[0]):] =[0,0,0,0,0]
            else:
                water_phase = find_water_phase(Planet['pressure'][i],Planet['temperature'][i])
                if water_phase == 'Water':

                    Phases[i][len(Mantle_phases[0]):] =[float('NaN'), 100, float('NaN'), float('NaN'),float('NaN')]
                elif water_phase == 'Ice_VII':

                    Phases[i][len(Mantle_phases[0]):] =[float('NaN'), float('NaN'), 100, float('NaN'),float('NaN')]
                elif water_phase == 'Ice_VI':
                    #Phases[i] =np.resize(Phases[i], len(Phases[i]) + 4)
                    Phases[i][len(Mantle_phases[0]):] = [float('NaN'), float('NaN'), float('NaN'), 100,float('NaN')]
                else:
                    Phases[i][len(Mantle_phases[0]):] = [float('NaN'), float('NaN'), float('NaN'), float('NaN'),100]

    else:
        if i < num_core_layers - 1:

            Phases[i][len(Mantle_phases[0]):] = [100]
        elif i < (num_mantle_layers + num_core_layers):
            # Phases[i] =np.resize(Phases[i],len(Phases[i])+4)

            Phases[i][len(Mantle_phases[0]):] = [0]
    return Phases



def get_alpha(Planet,core_wt_per,grids,layers):

    num_mantle_layers, num_core_layers, number_h2o_layers = layers

    total_layers = num_core_layers+num_mantle_layers+number_h2o_layers
    mantle_pressures = Planet['pressure'][num_core_layers:]
    mantle_temperatures = Planet['temperature'][num_core_layers:]
    core_pressures = Planet['pressure']
    core_temperatures = Planet['temperature']

    P_points_UM = []
    T_points_UM = []
    P_points_LM = []
    T_points_LM = []

    for i in range(len(mantle_pressures)):
        if mantle_pressures[i]<=1250000:
            P_points_UM.append(mantle_pressures[i])
            T_points_UM.append(mantle_temperatures[i])
        else:
            P_points_LM.append(mantle_pressures[i])
            T_points_LM.append(mantle_temperatures[i])

    alpha_UM = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']), grids[0]['alpha'],
                         (P_points_UM, T_points_UM), method='linear')

    alpha_LM = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']), grids[1]['alpha'],
                         (P_points_LM, T_points_LM), method='linear')

    alpha_mantle = np.concatenate((alpha_LM,alpha_UM),axis=0)

    #Core_speeds = minphys.get_core_speeds(core_pressures,core_temperatures,core_wt_per)
    #pdb.set_trace()
    alpha = []

    for i in range(total_layers):

        if i < num_mantle_layers + num_core_layers and i > num_core_layers:

            alpha.append(alpha_mantle[i - num_core_layers])
        else:
            alpha.append(0.)

    return alpha


def get_speeds(Planet,core_wt_per,grids,layers):
    num_mantle_layers, num_core_layers, number_h2o_layers = layers

    total_layers = num_core_layers+num_mantle_layers+number_h2o_layers
    mantle_pressures = Planet['pressure'][num_core_layers:]
    mantle_temperatures = Planet['temperature'][num_core_layers:]
    core_pressures = Planet['pressure']
    core_temperatures = Planet['temperature']

    P_points_UM = []
    T_points_UM = []
    P_points_LM = []
    T_points_LM = []

    for i in range(len(mantle_pressures)):
        if mantle_pressures[i]<=1250000:
            P_points_UM.append(mantle_pressures[i])
            T_points_UM.append(mantle_temperatures[i])
        else:
            P_points_LM.append(mantle_pressures[i])
            T_points_LM.append(mantle_temperatures[i])

    Mantle_speeds_UM = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']), grids[0]['speeds'],
                         (P_points_UM, T_points_UM), method='linear')

    Mantle_speeds_LM = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']), grids[1]['speeds'],
                         (P_points_LM, T_points_LM), method='linear')

    Mantle_speeds = np.concatenate((Mantle_speeds_LM,Mantle_speeds_UM),axis=0)

    Core_speeds = minphys.get_core_speeds(core_pressures,core_temperatures,core_wt_per)

    Vphi = []
    Vp = []
    Vs = []
    for i in range(total_layers):
        if i < num_core_layers:
            Vphi.append(Core_speeds[0][i]/1000.)
            Vp.append(Core_speeds[1][i]/1000.)
            Vs.append(Core_speeds[2][i]/1000.)

        elif i < num_mantle_layers + num_core_layers:
            Vphi.append(Mantle_speeds[i-num_core_layers][0])
            Vp.append(Mantle_speeds[i-num_core_layers][1])
            Vs.append(Mantle_speeds[i-num_core_layers][2])
        else:
            Vphi.append(0.)
            Vp.append(0.)
            Vs.append(0.)

    return Vphi,Vp,Vs

def write(Planet,filename):
    output = []
    for i in range(len(Planet['pressure'])):
        line_item = [(Planet['radius'][-1]-Planet['radius'][i])/1000.,Planet['radius'][i]/1000.,
                     Planet['density'][i]/1000.,Planet['pressure'][i]/10000.,Planet['temperature'][i],
                     Planet['Vphi'][i],Planet['Vp'][i],Planet['Vs'][i]]
        for j in range(len(Planet['phases'][i])):
            line_item.append(Planet['phases'][i][j])

        output.append(line_item)
    line_name = []
    line_name.append('Depth')
    line_name.append('Radius')
    line_name.append('Density')
    line_name.append('Pressure')
    line_name.append('Temperature')
    line_name.append('Vphi')
    line_name.append('Vp')
    line_name.append('Vs')
    for i in Planet['phase_names']:
        line_name.append(str(i))

    string_element = '  '.join(line_name)
    np.savetxt(filename+'.dat', output, '%.5f', "\t", newline='\n',
                header=string_element, footer='', comments='# ')

    print
    print "file written to:", filename+'.dat'
    print
    return 0

def find_water_phase(Pressure, Temperature):
    Pressure = Pressure/10000. #convert to GPa

    if Pressure > 20.3746:
        # ice VII
        phase = 'Ice_VII'

    elif Temperature < 355 and Temperature > 273.31:
        phase = iceVI(Temperature,Pressure)

    elif Temperature <= 715 and Temperature >= 355:
        # liquid or iceVII, run routine to check
        phase = iceVII(Temperature, Pressure)

    elif Pressure > 223.276 and Temperature < 355 and Temperature > 250:
        # iceVII or liq??
        # for bme3
        phase = iceVII_2(Temperature, Pressure)

    elif Pressure < 223.276 and Temperature < 355 and Temperature >= 273:
        # liquid water
        phase = 'Water'
    elif Pressure < 223.276 and Temperature < 273 and Temperature > 251:
        # liquid or iceIh, run routine to check
        phase = iceIh(Temperature, Pressure)

    elif Temperature <= 251 and Temperature > 0:
        # ice ih or ice VII
        phase = ih_or_vii(Pressure, Temperature)

    else:
        print "\n\n Outside Water phase diagram, need to quit \n"
        print "T = %r\tP = %r GPa" % (Temperature, Pressure)
        sys.exit()

    return phase

def iceVI(Temperature,Pressure):
    Theta = Temperature/273.31
    P_test = 632.4e6 * (1.-1.07476*(1.-pow(Theta,4.6)))
    if (Pressure*1e9) > P_test:
        Theta_VII = Temperature/355.
        P_test_VII = (1./1000.)*(2210.+534.2*((pow(Theta_VII,5.22)-1)))
        if Pressure > P_test_VII:
            return 'Ice_VII'
        else:
            return 'Ice_VI'
    else:
        return 'Water'

def iceVII(Temperature,Pressure):
    a = 1.73683
    b = 0.0544606
    c = 0.806106e-7
    Theta = Temperature/355.
    lnP = 2216.e6*np.exp(a*(1.-(1./Theta)) - b*(1.-pow(Theta,5.)) + c*(1.-pow(Theta,22.)))


    if (1.e9*Pressure)>lnP:
        return 'Ice_VII'

    else:
        return 'Water'

def iceVII_2(Temperature,Pressure):
    Tin = Temperature-250.
    Pm = 0.02186*Tin + 5.40841
    Pm = np.exp(Pm)*1.e6

    if (1.e9*Pressure) > Pm:
        print "in here"
        return 'Ice_VII'

    else:
        return 'Water'

def iceIH(Temperature,Pressure):
    a1 = 0.119539337e7
    a2 = 0.808183159e5
    a3 = 0.333826860e4
    b1 = 0.300000e1
    b2 = 0.257500e2
    b3 = 0.103750e3
    O  = Temperature/273.16
    pi = 1 + a1*(1.-O**b1) + a2*(1.-O**b2) + a3*(1.-O**b3)
    Pm = pi* 611.657

    if P > Pm:
        return 'Ice_Ih'
    else:
        return 'Water'

def Ih_or_VII(Temperature,Pressure):
    P1h = 176.0 + 0.918 * (Temperature - 198.5)

    P1h = P1h * 1e6

    if (1.e9*Pressure) > P1h:
        return 'Ice_VII'
    else:
        return 'Ice_Ih'




def update_radius(Planet, layers, rho_old):
    REarth = 6.371e6
    mass    = Planet.get('mass') #cumulative mass
    dmass   = Planet.get('dmass') #mass in each layer
    rad     = Planet.get('radius')
    grav    = Planet.get('gravity')
    rho_new = Planet.get('density')

    num_mantle_layers, num_core_layers, number_h2o_layers = layers
    n_tot = num_mantle_layers + num_core_layers + number_h2o_layers

    #density of each layer after dampening mixture
    rho   = np.zeros(n_tot)
    vol   = np.zeros(n_tot)

    mix   = 0.5
    delta = 0.0

    for i in range(0, n_tot):

        deltaNew = np.abs((rho_new[i]/rho_old[i]) - 1.0)
        if (deltaNew > delta):
            delta = deltaNew

        # the new density of this shell is a mix of new and old for dampening purposes
        rho[i] = rho_new[i] * mix + rho_old[i] * (1.0 - mix)

        # update radius and volume and then temperature
        if i > 0:
            vol[i] = dmass[i] / rho[i]
            rad[i] = (0.75 * (dmass[i] / rho[i] / np.pi) + (rad[i - 1] ** 3)) ** (1.0 / 3)
            #print 'r[%r] = %.3f' %(i,rad[i]/REarth)

    return(rad, delta)


def find_Planet_radius(radius_planet, core_mass_frac, structure_params, compositional_params, grids, Core_wt_per, layers):

    import planet
    import minphys

    def calc_CRF(value, args):
        radius_planet        = args[0]
        structure_params     = args[1]
        compositional_params = args[2]
        num_core_layers      = args[3][1]
        grids       = args[4]
        Core_wt_per = args[5]
        CMF_to_fit  = args[6]


        structure_params[6] = value
        Planet      = planet.initialize_by_radius(*[radius_planet, structure_params, compositional_params, layers])
        Planet      = planet.compress(*[Planet, grids, Core_wt_per, structure_params, layers])
        planet_mass = minphys.get_mass(Planet,layers)

        CMF = planet_mass[num_core_layers]/planet_mass[-1]
        print "Diff in Core Mass Fraction = ", '%.3e' % (CMF_to_fit - CMF)
        return (CMF_to_fit - CMF)

    def calc_CRF_WRF(values, *args):
        radius_planet        = args[0]
        structure_params     = args[1]
        compositional_params = args[2]
        num_core_layers      = args[3][1]
        num_mantle_layers    = args[3][0]
        num_water_laters     = args[3][2]
        grids       = args[4]
        Core_wt_per = args[5]
        CMF_to_fit  = args[6]
        WMF_to_fit  = args[7]


        structure_params[6] = values[0]
        structure_params[8] = values[1]
        #print values

        if (values[0]+values[1]) >= 1:
            if values[0] > values[1]:
                return (-200.,0)
            if values[1] > values[0]:
                return (0,-200.)

        Planet = planet.initialize_by_radius(*[radius_planet, structure_params, compositional_params, layers])
        Planet = planet.compress(*[Planet, grids, Core_wt_per, structure_params, layers])


        planet_mass      = minphys.get_mass(Planet,layers)
        core_mass        = planet_mass[num_core_layers]
        terrestrial_mass = planet_mass[num_core_layers+num_mantle_layers]
        water_mass       = planet_mass[-1] - terrestrial_mass

        CMF = core_mass/terrestrial_mass
        WMF = water_mass/planet_mass[-1]
        print '\n<========================================================================>'
        print "Diff in Core Mass percent = ", '%.3f' % (100.*CMF_to_fit - 100.*CMF)
        print "Diff in Water Mass percent = ", '%.3f' % (100.*WMF_to_fit - 100.*WMF)
        print '\ncalculated core mass frac {}'.format(CMF)
        print 'input core mass frac {}'.format(CMF_to_fit)
        print '\ncalculated h2o mass fraction {}'.format(WMF)
        print 'input h2o mass fraction {}'.format(WMF_to_fit)
        print '<========================================================================>'

        return (100.*CMF_to_fit - 100.*CMF,100.*WMF_to_fit - 100.*WMF)

    #if user enters water mass fraction and water layers
    if compositional_params[0] > 0:
        from scipy.optimize import root
        if layers[2] <10:
            print "***Build error: too few water layers for H2Owt% = {0:6f}***".format(compositional_params[0])
            print "Solution: changing to 100 water layers"
            layers[2] = 100

        water_mass_frac = compositional_params[0]
        args     = (radius_planet, structure_params, compositional_params, layers, grids, Core_wt_per, core_mass_frac,water_mass_frac)
        solution = root(calc_CRF_WRF, [.5,water_mass_frac],args=args,tol=1.e-4,method='anderson')
        structure_params[6], structure_params[8] = solution.x
        if structure_params[6]> 1.0:
            structure_params[6] = 0.5
        print 'core rad frac {}'.format(structure_params[6])
        print 'h2o rad fraction {}'.format(structure_params[8])

        Planet = planet.initialize_by_radius(*[radius_planet, structure_params, compositional_params, layers])
        Planet = planet.compress(*[Planet, grids, Core_wt_per, structure_params, layers])
        return Planet

    else:
        from scipy.optimize import brentq
        #make sure there are no water layers. User has entered H2Owt% = 0
        #if user entered >0 water layers or >0 water radius fraction, zero them out
        if layers[2] > 0 or structure_params[8] != 0 :
            print "***Build error: excess in water layers for H2Owt% = 0 wt%***"
            print "Solution: removing water layer"
            structure_params[8] = 0
            layers[2] = 0

        args = [radius_planet, structure_params, compositional_params, layers,grids,Core_wt_per,core_mass_frac]
        structure_params[6] = brentq(calc_CRF,.30,.85,args=args,xtol=1e-4)

        Planet = planet.initialize_by_radius(*[radius_planet, structure_params, compositional_params, layers])
        Planet = planet.compress(*[Planet, grids, Core_wt_per, structure_params, layers])


        return Planet


def R_of_M(mass_planet, core_mass_frac, structure_params, compositional_params, grids, Core_wt_per, layers):
    import planet
    import minphys
    #1. initialize the plant by mass
    Planet = planet.initialize_by_mass(*[mass_planet, structure_params, compositional_params, layers, core_mass_frac])

    #DEBUG
    #print 'Mass after initialization = {}'.format(Planet.get('mass')[-1]/5.972e24)
    #print 'Surface pressure = {}'.format(Planet.get('pressure')[-1])
    #print mass_planet
    #sys.exit()


    #2. run routine to find density as a function of P, T from perplex files
    # iterate until the max change in density of any given layer is 1e-6

    Planet = planet.compress_fixed_mass(*[Planet, grids, Core_wt_per, structure_params, layers])




    print 'Mass after compression = {}'.format(Planet.get('mass')[-1]/5.972e24)

    #Planet['phases'] = functions.get_phases(Planet, grids, layers)

    #Planet['phase_names'] = names

    #Planet['Vphi'], Planet['Vp'], Planet['Vs'] = functions.get_speeds(Planet, Core_wt_per, grids, layers)


    return Planet





