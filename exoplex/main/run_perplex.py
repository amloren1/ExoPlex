import os
import sys
import pexpect as pe
import pdb

PerPlex_path = os.path.dirname(os.path.realpath(__file__))+"/PerPlex"
cur_path = os.path.dirname(os.path.realpath(__file__))

# hack to allow scripts to be placed in subdirectories next to ExoPlex:
if not os.path.exists('main') and os.path.exists('../main'):
    
    sys.path.insert(1, os.path.abspath('..'))
    


if not os.path.exists(PerPlex_path+'/build') or not os.path.exists(PerPlex_path+'/vertex') \
    or not os.path.exists(PerPlex_path+'/werami'):

    import urllib2
    import tarfile

    print '\n\n** Downloading Perple_x. Please wait**'
    print 'url: {}'.format('http://www.perplex.ethz.ch/ExoPlex/')

    perplex_linux = 'http://www.perplex.ethz.ch/ExoPlex/Perple_X_6.8.1_Linux_64_gfortran_edited.tar.gz'
    perplex_mac   = 'http://www.perplex.ethz.ch/ExoPlex/Perple_X_6.8.1_OSX_10.6+_Intel_MC_Mar_6_2018_edited.zip'

    platform = sys.platform
    
    #check platform to download correct perple_x version
    if platform != 'linux2':
        perplex_link = perplex_mac
    else:
        perplex_link = perplex_linux
        
    try:
        filename = 'perplex.tar.gz'
        response = urllib2.urlopen(perplex_link)
        
        with open(filename, 'wb') as f: f.write(response.read())
        
        tar = tarfile.open(filename, 'r:gz')
        
        tar.extract('build', PerPlex_path)
        tar.extract('vertex', PerPlex_path)
        tar.extract('werami', PerPlex_path)
        print '\nSuccessfully downloaded Perple_x. Stored in {}\n\n'.format(PerPlex_path)

    except:
        print '\nUnable to download Perple_x. Make sure you are connected to the internet'
        
        
    


def download_perplex():
    import urllib2
    import tarfile
    
    print '\n\n** Downloading Perple_x. Please wait**'
    print 'url: {}'.format('http://www.perplex.ethz.ch/perplex/exoplex/')
    
    perplex_linux = 'http://www.perplex.ethz.ch/perplex/exoplex/linux/Perple_X_6.8.1_Linux_64_gfortran.tar.gz'
    perplex_mac   = 'http://www.perplex.ethz.ch/perplex/exoplex/OSX/Perple_X_6.8.1_OSX_10.6+_Intel_MC_Mar_6_2018.zip'
    
    try:
        filename = 'perplex.tar.gz'
        response = urllib2.urlopen(perplex_linux)
        
        with open(filename, 'wb') as f: f.write(response.read())
        
        tar = tarfile.open(filename, 'r:gz')
        
        tar.extract('build', PerPlex_path)
        tar.extract('vertex', PerPlex_path)
        tar.extract('werami', PerPlex_path)
        
    except:
        print 'Unable to download Perple_x. Make sure you are connected to the internet'
        return
    print 'Successfully downloaded Perple_x. Stored in {}'.format(PerPlex_path)
    
    return
    


def run_perplex(*args):
    

    Mantle_wt_per = args[0]


    solutionFileNameMan = args[1]

    Pressure_range_mantle = args[2][0]
    Temperature_range_mantle = args[2][1]
    resolution = args[2][2]

    filename = args[3]
    UMLM = args[4]

    plxMan = str(Mantle_wt_per.get('MgO')) + ' ' + str(Mantle_wt_per.get('SiO2')) + ' ' \
             + str(Mantle_wt_per.get('FeO')) + ' ' + str(Mantle_wt_per.get('CaO')) \
             + ' ' + str(Mantle_wt_per.get('Al2O3'))+ ' ' + str(0.) #last value included for Na



    filename ='Solutions/'+solutionFileNameMan

    if os.path.isfile(filename+'_UM.tab') and UMLM == True:
        print '\nThe Upper mantle .tab already exists, please wait briefly for solution:'
        print  filename+'_UM.tab\n'
        return filename

    if os.path.isfile(filename+'_LM.tab') and UMLM == False:
        print '\nThe Lower mantle .tab already exists, please wait briefly for solution'
        print  filename+'_LM.tab\n'
        return filename

    else:
        if UMLM == True:
            print 'Making upper mantle PerPlex phase file. \n This will be stored in: '+ filename+'_UM.tab'
            description = '_UM'
        else:
            print 'Making lower mantle PerPlex phase file. \n This will be stored in: '+ filename+'_LM.tab'
            description = '_LM'

        #we need to shorten the file name for PerPlex to accept it
        solutionFileNameMan_short = list(filename)
        solutionFileNameMan_short[0:30] = []

        solutionFileNameMan = "".join(solutionFileNameMan_short)+description


    # define perplex inputs in terms of components, this is for legibility

    component1 = 'MGO'
    component2 = 'SIO2'
    component3 = 'FEO'
    component4 = 'CAO'
    component5 = 'AL2O3'
    component6 = 'NA2O'

    p = pe.spawn(PerPlex_path+"/./build")


    p.sendline(solutionFileNameMan)
    p.sendline(PerPlex_path+'/stx11ver.dat')
    p.sendline(PerPlex_path+'/perplex_option.dat')
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
    p.sendline(Temperature_range_mantle)
    # Enter minimum and maximum values, respectively, for: P(bar)
    # P(Pa) = P(bar)*100000
    p.sendline(Pressure_range_mantle)
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
    p.sendline(PerPlex_path+'/solution_model.dat')
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
    p.sendline('Gt') #kitchen sink
    p.sendline('Ppv')
    p.sendline('CF')
    p.sendline('')
    # Enter calculation title:
    p.sendline(solutionFileNameMan + 'calc')

    p.logfile = open('build.log','wb')
    p.read()
    p.wait()

    print "Done with Build, moving on to Vertex"

    # Spawn Vertex ----------------#
    # Enter the project name (the name assigned in BUILD) [default = my_project]:
    p = pe.spawn(PerPlex_path+"/./vertex",timeout=2400)

    p.sendline(solutionFileNameMan)

    #p.expect('$$$$',timeout=None)
    p.logfile = open('vertex.log', 'wb')
    p.read()
    p.wait()

    
    print 'Finished with Vertex, beginning Werami'

    try:

        p = pe.spawn(PerPlex_path+"/./werami",timeout=None)


        p.sendline(solutionFileNameMan)
        # select 2D grid
        p.sendline('2')
        # Below, select parameters density, alpha, cp.
        # Ns for no calculating individual phase properties
        p.sendline('2')
        p.sendline('N')
        p.sendline('12')
        p.sendline('N')
        p.sendline('13')
        p.sendline('N')
        p.sendline('14')
        p.sendline('N')
        p.sendline('4')
        p.sendline('N')
        p.sendline('19')
        p.sendline('N')
        ####### the next lines will pass requests to perplex to print phases and their proportions into the .tab file
        phases = ['C2/c', 'Wus', 'Pv', 'an', 'Sp', 'O', 'Wad', \
        'Ring', 'Opx', 'Cpx', 'Aki', 'Gt', 'Ppv', 'CF', 'st', \
        'q', 'ca-pv', 'cfs', 'coe', 'ky', 'seif' ]
        
        p.sendline('7')
        for phase in phases:
            p.sendline(phase)
            ii = p.expect(['try again:','proportions keyword'])
            if ii == 0:
                #print phase
                continue
                #sys.exit()
            elif phase != 'seif' and ii == 1:
                p.sendline('7')
                continue
            else:
                continue
       
        # exit parameter choosing

        p.sendline('0')
        # Change default variable range (y/n)?
        p.sendline('N')

        # Enter number of nodes in the T(K)     and P(bar)   directions:
        p.sendline(resolution)
        
        p.logfile = open('werami.log','wb')
        p.expect('EXIT', timeout=None)
        p.terminate()
        print "Done with PerPlex"

        if UMLM == True:
            os.rename(solutionFileNameMan+'_1.tab', filename+'_UM.tab')
        else:
            os.rename(solutionFileNameMan + '_1.tab', filename + '_LM.tab')

        successful = True
    except:

        successful = False
        print 'perplex broke at werami. The details are stored in the ERROR_ files'
        os.rename('build.log', 'ERROR_'+solutionFileNameMan+'_build.log')
        os.rename('vertex.log', 'ERROR_'+solutionFileNameMan+'_vertex.log')
        os.rename('werami.log', 'ERROR_'+solutionFileNameMan+'_werami.log')



    os.remove(solutionFileNameMan+'.arf')
    os.remove(solutionFileNameMan+'.blk')
    os.remove(solutionFileNameMan+'.dat')
    os.remove(solutionFileNameMan+'.plt')
    os.remove(solutionFileNameMan+'.tof')
    #os.remove(solutionFileNameMan+'_VERTEX_options.txt')
    #os.remove(solutionFileNameMan+'_WERAMI_options.txt')
    #os.remove(solutionFileNameMan+'_auto_refine.txt')

    return filename


    
    
    
'''

 p.sendline('q')  # 13
        p.sendline('7')
        p.sendline('ca-pv')  # 14
        p.sendline('7')
        p.sendline('cfs')  # 15
        p.sendline('7')
        p.sendline('coe')  # 16
        p.sendline('7')
        p.sendline('ky')  # 17
        p.sendline('7')
        p.sendline('seif')  # 18

        # 21 species, in all for Fe-Si-Mg-O regime
        p.sendline('7')
        p.sendline('C2/c')  # 0
        p.sendline('7')
        p.sendline('Wus')  # 1
        p.sendline('7')
        p.sendline('Pv')  # 2
        p.sendline('7')
        p.sendline('an')  # 3
        p.sendline('7')
        p.sendline('Sp')  #4
        p.sendline('7')
        p.sendline('O')  # 4
        p.sendline('7')
        p.sendline('Wad')  # 5
        p.sendline('7')
        p.sendline('Ring')  # 6  
        p.sendline('7')
        p.sendline('Opx')  # 7
        p.sendline('7')
        p.sendline('Cpx')  # 8
        p.sendline('7')
        p.sendline('Aki')  # 9
        p.sendline('7')
        p.sendline('Gt')  # 10 gt_maj
        p.sendline('7')
        p.sendline('Ppv')  # 11
        p.sendline('7')
        p.sendline('CF')   # 12
        p.sendline('7')
        p.sendline('st')  # 12
        p.sendline('7')
'''
