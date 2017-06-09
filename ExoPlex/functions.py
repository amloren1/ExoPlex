import numpy as np
import sys

def get_percents(*args):
    FeMg = args[1]
    SiMg = args[2]
    CaMg = args[3]
    AlMg = args[4]
    mol_frac_Fe_mantle = args[5]
    wt_frac_Si_core = args[6]
    wt_frac_O_core = args[7]
    wt_frac_S_core = args[8]

    MgSi = SiMg**(-1)
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

    #Sum of all masses = 1 g
    b = np.array([0., 0. , 0. , 0. ,  0. , 0. , 0 , 0.,0., 100.])

############################################################################################

    #Comment here
    #Double checks
    #Return only moles of elements, no molecules
    A = np.array([ [0., 0., 0. , 0. , -1. ,-1, -2 , 1 , -1, -1.5] ,
            [1 , -FeSi , 0. , 0. , 1 , 0. , -FeSi , 0.,0.,0.] ,
           [0, -MgSi, 0, 0 ,0. ,1 ,-MgSi ,0.,0.,0.] ,
            [0, -CaSi, 0, 0 ,0. ,0. ,-CaSi ,0., 1 ,0.] ,
        [0, -AlSi, 0, 0 ,0. ,0. ,-AlSi ,0.,0., 1] ,
          [mol_frac_Fe_mantle , 0 , 0 ,0, (mol_frac_Fe_mantle-1) ,  0 , 0 , 0 ,0., 0] ,
          [wt_frac_Si_core*mFe , (wt_frac_Si_core-1)*mSi , wt_frac_Si_core*mO , wt_frac_Si_core*mS \
                       , 0. , 0 , 0. , 0 , 0., 0.] ,
          [wt_frac_O_core*mFe , (wt_frac_O_core)*mSi , (wt_frac_O_core-1)*mO , wt_frac_O_core*mS
                       , 0. , 0 , 0. , 0 , 0., 0.]  ,
          [wt_frac_S_core*mFe , (wt_frac_S_core)*mSi , (wt_frac_S_core)*mO , (wt_frac_S_core-1)*mS
                       , 0. , 0 , 0. , 0 , 0., 0.] ,
          [mFe , mSi , mO, mS ,  mFe , mMg , mSi ,mO, mCa, mAl]   ])


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

    mass_of_Core = (mFe*(Core_moles[0])+mSi*Core_moles[1]\
                    +(mO)*Core_moles[2]+mS*Core_moles[3])

    mass_of_Mantle = (mFe*Mantle_moles[0]+mMg*Mantle_moles[1]\
                   +mSi*Mantle_moles[2]+mO*Mantle_moles[3]\
                   +mCa*Mantle_moles[4]+mAl*Mantle_moles[5])

    Mtot= mass_of_Core+mass_of_Mantle #in g

    #Weight percents of mantle oxides
    #Weight percents assuming FeO, MgO, SiO2, CaO, Al2O3

    FeO_mant_wt = Mantle_moles[0]*(mFe+mO)/mass_of_Mantle
    MgO_mant_wt = Mantle_moles[1]*(mMg+mO)/mass_of_Mantle
    SiO2_mant_wt = Mantle_moles[2]*(mSi+2.*mO)/mass_of_Mantle
    CaO_mant_wt = Mantle_moles[4]*(mCa+mO)/mass_of_Mantle
    Al2O3_mant_wt = (Mantle_moles[5]/2.)*(2.*mAl+3.*mO)/mass_of_Mantle



    #Throw exception, not if statement
    #make inequality not, absolute if. Use machine precision
    if (FeO_mant_wt+MgO_mant_wt+SiO2_mant_wt+CaO_mant_wt+Al2O3_mant_wt) != 1.:
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

    Mantle_wt_per = {'FeO': FeO_mant_wt, 'SiO2': SiO2_mant_wt, 'MgO': MgO_mant_wt, 'CaO': CaO_mant_wt,'Al2O3':Al2O3_mant_wt}

    #decimal fraction of materials in CORE by mass, these are perplex inputs (hence the rounding)
    # this is the wt% of Fe and the Fe in FeSi *** so total Fe in core

    Fe_core_wt   = (Core_moles[0])*(mFe)/mass_of_Core
    Si_core_wt   = Core_moles[1]*(mSi)/mass_of_Core
    O_core_wt   = Core_moles[2]*(mO)/mass_of_Core
    S_core_wt   = Core_moles[3]*(mS)/mass_of_Core

    #Throw exception, not if statement
    #make inequality not, absolute if. Use machine precision
    if (S_core_wt+O_core_wt+Si_core_wt+Fe_core_wt) != 1.:
        print '\n\n damn, it broke'
        sys.exit()

    Fe_core_wt = abs(round(Fe_core_wt*100.,8))
    Si_core_wt  = abs(round(Si_core_wt*100.,8))
    O_core_wt = abs(round(O_core_wt*100.,8))
    S_core_wt  = abs(round(S_core_wt*100.,8))

    Core_wt_per = {'Fe':Fe_core_wt,'Si':Si_core_wt,'O':O_core_wt,'S':S_core_wt}
    Core_mol_per ={'Fe':Core_moles[0]/tot_moles_core,'Si':Core_moles[1]/tot_moles_core,\
                  'O':Core_moles[2]/tot_moles_core,'S':Core_moles[3]/tot_moles_core}


    return(Core_wt_per,Mantle_wt_per,Core_mol_per)


def verbosity():

    # MgO,SiO2,FeO, CaO, Al2O3 wt%


    solutionFileNameMan    = 'SiMg_FeMg_CaMg_AlMg_XFeO_fSic_fOc_fSc' + solfileparamsString+'_MANTLE'
#    datafileName         = 'Mass_MgSi_FeSi_XFeO_fSic_fOc_' +repr(round(PlanetMass/MEarth,3))+solfileparamsString

    print 'Mantle solution file: \n'
    print solutionFileNameMan

    #sys.exit()
    #for creating data files in perplex, plcCor is an entry into the build code
    plxCor  = repr(Fecwt)+' '+repr(O2cwt)+' '+repr(Sicwt) + ' ' + repr(S2cwt)
    solfileparamsString0 = '_'+repr(round(fSic,3)) +'_'+repr(round(fOc,3)) + '_'+ repr(round(fSc,3))
    solfileparamsString  = solfileparamsString0.replace('.',',')
#    datafileName         = 'Mass_MgSi_FeSi_XFeO_fSic_fOc_' +repr(normMass)+solfileparamsString
    solutionFileNameCor     = 'fSic_fOc_fSc' + solfileparamsString+'_CORE'

     #print some stuff here
    if verbose:
        print 'Mantle composition (wt%%): \n'
        print 'MgO = %.4f \tSiO2 = %.4f \tFeO = %.4f \tCaO = %.4f \tAl2O3 = %.4f\n' %                        \
                (np.abs(MgOmwt),np.abs(SiO2mwt),np.abs(FeOmwt),np.abs(CaOmwt), np.abs(Al2O3mwt) )

        print 'Core composition (wt%%): \n'
        print 'O2 = %.4f \tSi = %.4f \tS2 = %.4f \tFe = %.4f (Fe alone & Fe in FeSi)\n' %          \
             (np.abs(O2cwt),np.abs(Sicwt),np.abs(S2cwt),np.abs(Fecwt))

        print 'wt%% of Core: %r \t wt%% of Mantle: %r' % (round(masfCor*100,5), round(masfMan*100,5))
        print '*---------------------------------------------------*\n'



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
    # x = [nFec,nMgc,nSic,nOc, nSc | ,nFem,nMgm,nSim,nOm]
    if verbose:
        print 'Calculated Values'
        print 'Mg/Si = %.4f' % (x[5]/(x[1]+x[6]))
        print 'Fe/Si = %.4f' % ((x[0]+x[4])/(x[1]+x[6]))
        print 'Ca/Si = %.4f' % ((x[8])/(x[1]+x[6]))
        print 'Al/Si = %.4f' % ((x[9])/(x[1]+x[6]))
        print 'XFeO = %.4f' % (x[4]/(x[0]+x[4]))
        print 'fSic = %.4f' % ((x[1]*(mSi))/masfCor)
        print 'fOc = %.4f' % ((x[2]*(mO))/masfCor)
        print 'fSc = %.4f' % ((x[3]*(mS))/masfCor)
        print 'Core wt%% of Si = %.4f %%' % (100*(x[1]*(mSi))/masfCor)
        print 'Core wt%% of O = %.4f %%' % (100*(x[2]*(mO))/masfCor)
        print 'Core wt%% of S = %.4f %%' % (100*(x[3]*(mO))/masfCor)
        print 'Total light element wt%% in core: %.4f %%' % (100*(x[2]*(mO)+x[1]*mSi+x[3]*mS)/masfCor)
        print 'Mantle Mg/Si = %.4f' % (x[5]/x[6])
        print '\nEntered values:'
        print 'Si/Mg = %r' % SiMg
        print 'Fe/Mg = %r' % FeMg
        print 'Ca/Mg = %r' % CaMg
        print 'Al/Mg = %r' % AlMg
        print 'XFeO = %.4f' % XFeO
        print 'fSic = %.4f' % (fSic)
        print 'fOc = %.4f' % fOc
        print 'fSc = %.4f' % fSc

    #sys.exit()
#  return values: masfCor,masfMan, MgOmwt,SiO2mwt,FeOmwt, Sicwt, Fecwt,O2cwt
    #must rememer: the masfMan and Cor are wt%s
    return (masfCor,masfMan,MgOmwt,SiO2mwt,FeOmwt,CaOmwt, Al2O3mwt, Sicwt,Fecwt,O2cwt,S2cwt,solutionFileNameCor,solutionFileNameMan, plxMan,plxCor)

