import os
import sys

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))

def run_perplex(Mantle_wt_per,SiMg,FeMg,CaMg,AlMg,mol_frac_Fe_mantle,wt_frac_Si_core,wt_frac_O_core,wt_frac_S_core,verbose):
    solutionFileNameMan = 'test'

    if os.path.isfile('../Solutions/'+solutionFileNameMan):
        print "true"
        if verbose:  # the verbose variable is
            print '\n\n***The mantle .tab already exists, please wait briefly for solution***\n'
            print 'Mantle File name: ' + pman + solutionFileNameMan
        return

    return 0
"""
    else:
        if verbose == 'y':
            print '\n\n***new mantle solution file will be generated:\n', solutionFileNameMan, ' ***'
            print 'it will be stored in: ', pman

            sman = list(solutionFileNameMan)
            sman[0:38] = []

            solutionFileNameMan = "".join(sman)
            # print pman+solfile
            # sys.exit()

    plxMan   = str(Mantle_wt_per.get('MgO')) + ' ' + str(Mantle_wt_per.get('SiO2')) + ' ' \
               +str(Mantle_wt_per.get('FeO'))+ ' ' +str(Mantle_wt_per.get('CaO')) \
                + ' ' +str(Mantle_wt_per.get('Al2O3'))




    solfileparamsString0 = '_'+str(round(SiMg,3))+'_'+str(round(FeMg,3))+'_'+str(round(CaMg,3))+'_'+ str(round(AlMg,3))\
                           +'_'+str(round(mol_frac_Fe_mantle,3)) +'_'+str(round(wt_frac_Si_core,3)) \
                           +'_'+str(round(wt_frac_O_core,3)) + '_'+ str(round(wt_frac_S_core,3))

    solfileparamsString  = solfileparamsString0.replace('.',',')




    ## define perplex inputs ##
    component1    = 'MGO'
    component2    = 'SIO2'
    component3    = 'FEO'
    component4    = 'CAO'
    component5    = 'AL2O3'


    return ()
    """