import os
import sys
import pexpect as pe

PerPlex_path = os.path.dirname(os.path.realpath(__file__))+"/PerPlex"

# hack to allow scripts to be placed in subdirectories next to ExoPlex:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))

def run_perplex(Mantle_wt_per,SiMg,FeMg,CaMg,AlMg,mol_frac_Fe_mantle,wt_frac_Si_core, prange, trange,verbose):


    plxMan = str(Mantle_wt_per.get('MgO')) + ' ' + str(Mantle_wt_per.get('SiO2')) + ' ' \
             + str(Mantle_wt_per.get('FeO')) + ' ' + str(Mantle_wt_per.get('CaO')) \
             + ' ' + str(Mantle_wt_per.get('Al2O3'))+ ' ' + str(0.) #last value included for Na


    solfileparamsString0 = '_' + str(round(SiMg, 3)) + '_' + str(round(FeMg, 3)) + '_' + str(
        round(CaMg, 3)) + '_' + str(round(AlMg, 3)) \
                           + '_' + str(round(mol_frac_Fe_mantle, 3)) + '_' + str(round(wt_frac_Si_core, 3))

    # changes periods to commas
    solfileparamsString = solfileparamsString0.replace('.', ',')

    solutionFileNameMan = 'SiMg_FeMg_CaMg_AlMg_XFeO_fSic' + solfileparamsString + '_MANTLE'

    if os.path.isfile('../Solutions/'+solutionFileNameMan):
        if verbose:  # the verbose variable is
            print '\n\n***The mantle .tab already exists, please wait briefly for solution***\n'
            print 'Mantle File name: '+ solutionFileNameMan
        return

    else:
        print '\n\n***new mantle solution file will be generated:\n', solutionFileNameMan, ' ***'
        print 'it will be stored in: ../Solutions/'+ solutionFileNameMan
        #we need to shorten the file name for PerPlex to accept it
        solutionFileNameMan_short = list(solutionFileNameMan)
        solutionFileNameMan_short[0:30] = []

        solutionFileNameMan = "".join(solutionFileNameMan_short)
    # define perplex inputs in terms of components, this is for legibility

    component1 = 'MGO'
    component2 = 'SIO2'
    component3 = 'FEO'
    component4 = 'CAO'
    component5 = 'AL2O3'
    component6 = 'NA2O'

    p = pe.spawn(PerPlex_path+"/./build")


    p.sendline(solutionFileNameMan)
    p.sendline('../ExoPlex/PerPlex/stx11ver.dat')
    p.sendline('../ExoPlex/PerPlex/perplex_options.dat')
    # Transform them (Y/N)?
    p.sendline('N')
    # Calculations with saturated components (Y/N)?
    p.sendline('N')
    # Use chemical potentials, activities or fugacities as independent variables (Y/N)?
    p.sendline('N')
    # Select thermodynamic components from the set:
    p.sendline(component1)  # MGO
    p.sendline(component2)  # SIO2
    p.sendline(component3)  # FEO
    p.sendline(component4)  # CAO
    p.sendline(component5)  # AL2O3
    p.sendline(component6)  # NA2O
    p.sendline('')
    # Specify computational mode:

    p.sendline('2')

    # Make one dependent on the other, e.g., as along a geothermal gradient (y/n)?
    p.sendline('N')

    # Select x-axis variable:
    p.sendline('2')
    # Enter minimum and maximum values, respectively, for: T(K)
    p.sendline(trange)
    # Enter minimum and maximum values, respectively, for: P(bar)
    # P(Pa) = P(bar)*100000
    p.sendline(prange)
    # Specify component amounts by weight (Y/N)?
    p.sendline('Y')
    # Enter weight amounts of the components:
    # MGO SIO2 FEO CAO AL2O3
    # for the bulk composition of interest:
    # NOTE*** This is a wt%
    p.sendline(plxMan)
    # Output a print file (Y/N)?
    p.sendline('N')
    # Exclude pure and/or endmember phases (Y/N)?
    p.sendline('N')

    # Include solution models (Y/N)?
    p.sendline('Y')
    p.sendline('../ExoPlex/PerPlex/stx11_solution_model.dat')
    p.sendline('C2/c') #C2C Phase of clinopyroxene
    p.sendline('Wus')
    p.sendline('Pv')
    p.sendline('Pl')
    p.sendline('Sp')
    p.sendline('O')
    p.sendline('Wad')
    p.sendline('Ring')
    p.sendline('Opx')
    p.sendline('Cpx')
    p.sendline('Aki')
    p.sendline('Gt_maj') #kitchen sink
    p.sendline('Ppv')
    p.sendline('CF')
    p.sendline('')
    # Enter calculation title:
    p.sendline(solutionFileNameMan + 'calc')

    p.logfile = open('death.log','wb')
    p.read()
    p.wait()

    print "Done with Build, moving on to Vertex"

    # Spawn Vertex ----------------#
    # Enter the project name (the name assigned in BUILD) [default = my_project]:
    p = pe.spawn(PerPlex_path+"/./vertex",timeout=600)

    p.sendline(solutionFileNameMan)

    p.expect('$$$$',timeout=None)
    p.logfile = open('death.txt', 'wb')
    p.read()
    p.wait()
