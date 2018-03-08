import sys
import numpy as np
Earth_radius = 6.371e6
Earth_mass   = 5.972e24 # kg 
import minphys
import pdb
def initialize_by_radius(*args):
    radius_planet        = args[0]
    structural_params    = args[1]
    compositional_params = args[2]
    num_mantle_layers, num_core_layers, number_h2o_layers = args[3]

    wt_frac_water, FeMg, SiMg, CaMg, AlMg, mol_frac_Fe_mantle, wt_frac_Si_core, \
     wt_frac_O_core, wt_frac_S_core = compositional_params

    core_rad_frac         = structural_params[6]
    Mantle_potential_temp = structural_params[7]
    water_rad_frac        = structural_params[8]
    water_potential_temp  = structural_params[9]

    #def Setup(Mp, masMan, masCor, masH2O, nr, tCrstMan, dtManCor):
    # global Th,Tss
    #    global rad, rho, rho0, vol, dmass, mass, delta,\
    #        T, nic, noc, grav, P, rhon, alpha, cp,phase

    # use inputs to find actual mass of mantle and core

    # if there is a water layer, the imput temperature is lowered because that temperature is for the crustal layer
    # also 50 shells are used for the water layer hence the nh20 vaiable


    num_layers = num_core_layers+num_mantle_layers + number_h2o_layers

    # these grids are initialized in this function and are then a passed around through the other routines and edited each iteration
    radius_layers      = np.zeros(num_layers)
    density_layers     = np.zeros(num_layers)
    volume_layers      = np.zeros(num_layers)
    mass_layers        = np.zeros(num_layers)
    cumulative_mass    = np.zeros(num_layers)
    Temperature_layers = np.zeros(num_layers)
    # used for compressiofh2on funciton
    gravity_layers  = np.zeros(num_layers)
    Pressure_layers = np.zeros(num_layers)
    alpha           = np.zeros(num_layers)
    cp              = np.zeros(num_layers)
    Vphi            = np.zeros(num_layers)

    Vs = np.zeros(num_layers)
    Vp = np.zeros(num_layers)
    K  =  np.zeros(num_layers)

    planet_radius_guess    = radius_planet*Earth_radius
    water_thickness_guess  = water_rad_frac*planet_radius_guess
    core_thickness_guess   = core_rad_frac * (planet_radius_guess-water_thickness_guess)
    mantle_thickness_guess = planet_radius_guess - water_thickness_guess - core_thickness_guess



    for i in range(num_layers):
        if i < num_core_layers:
            radius_layers[i]      = ((float(i)/(num_core_layers-1.))*core_thickness_guess)
            density_layers[i]     = 8280.0
            Temperature_layers[i] = 4500.+ (float((3000.-4500.))/float(num_core_layers-1.))\
                                                            *float((i))

        elif i < (num_core_layers+num_mantle_layers):
            radius_layers[i]      = (core_thickness_guess+((float(i-num_core_layers)/(num_mantle_layers-1))*mantle_thickness_guess))
            density_layers[i]     = 4100.
            Temperature_layers[i] = 2700.

        elif number_h2o_layers>0:
            radius_layers[i]      = core_thickness_guess+mantle_thickness_guess+\
                                   ((float(i-num_core_layers-num_mantle_layers)/(number_h2o_layers-1.))*water_thickness_guess)
            density_layers[i]     = 1100.
            Temperature_layers[i] = 300.

    for i in range(num_layers):
        if i > num_core_layers+num_mantle_layers:
            Pressure_layers[i] = 1
        else:
            Pressure_layers[i] = (float((5000.-(300.*10000))/float(num_core_layers+num_mantle_layers))*float(i)
                                  + 300.*10000)


    keys = ['radius','density','temperature','gravity','pressure',\
            'alpha','cp','Vphi''Vp','Vs','K']


    return dict(zip(keys,[radius_layers, density_layers,Temperature_layers,gravity_layers, Pressure_layers,
                          alpha, cp,Vphi,Vp,Vs,K]))



def compress(*args):
    Planet            = args[0]
    grids             = args[1]
    Core_wt_per       = args[2]
    structural_params = args[3]
    layers            = args[4]

    n_iterations      = 1
    max_iterations    = 100
    n_min             = 2


    old_rho = [0  for i in range(len(Planet['density']))]
    converge = False
    print
    while n_iterations <= max_iterations and converge == False:
        print "iteration #",n_iterations


        #find density with current P, T gradients
        Planet['density'] = minphys.get_rho(Planet,grids,Core_wt_per,layers)

        for i in range(len(Planet['density'])):
            if np.isnan(Planet['density'][i]) == True:
                print "Density has a nan"
                print i, Planet['pressure'][i],Planet['temperature'][i]
                print
                sys.exit()
            elif np.isnan(Planet['gravity'][i]) == True:
                print "Gravity has a nan \nat shell\tPressure\tTemperature"
                print i, Planet['pressure'][i],Planet['temperature'][i]
                print
                sys.exit()

        #update gravity, pressure and temperature with new density
        Planet['gravity']     = minphys.get_gravity(Planet,layers)

        Planet['pressure']    = minphys.get_pressure(Planet,layers)

        #print Planet['pressure']
        Planet['temperature'] = minphys.get_temperature(Planet,grids,structural_params,layers)
        

        converge,old_rho = minphys.check_convergence(Planet['density'],old_rho)

        if n_iterations < n_min:
            converge = False
        n_iterations+=1

    return Planet
    
def compress_fixed_mass(*args):
    
    import functions as f
    import copy
    Planet            = args[0]
    grids             = args[1]
    Core_wt_per       = args[2]
    structural_params = args[3]
    layers            = args[4]

    n_iterations      = 1
    max_iterations    = 100
    n_min             = 2

    
    converge = False
    print
    while n_iterations <= max_iterations and converge == False:
        
        
        print "iteration #",n_iterations
        print 'radius = {}'.format(Planet.get('radius')[-1]/Earth_radius)
        
        #find density with current P, T gradients
        
        #copy the previous density before updating again
        plan = copy.deepcopy(Planet)
        old_rho = plan['density']
        import pdb
        #pdb.set_trace()
        #find new density
        Planet['density'] = minphys.get_rho(Planet,grids,Core_wt_per,layers)
        
        
        for i in range(len(Planet['density'])):
            if np.isnan(Planet['density'][i]) == True:
                print "Density has a nan"
                print i, Planet['pressure'][i],Planet['temperature'][i]
                print
                sys.exit()
            elif np.isnan(Planet['gravity'][i]) == True:
                print "Gravity has a nan \nat shell\tPressure\tTemperature"
                print i, Planet['pressure'][i],Planet['temperature'][i]
                print
                sys.exit()
        
        update = f.update_radius(Planet, layers, old_rho)

        Planet['radius']      = update[0]
        delta = update[1]
        
       
        #update gravity, pressure, and temperature using minphys,calc_ routines :)
        Planet['gravity']     = minphys.calc_gravity(Planet, layers)
        
        Planet['pressure']    = minphys.calc_pressure(Planet, layers)
        
        #pdb.set_trace()
        
        Planet['temperature'] = minphys.get_temperature(Planet,grids,structural_params,layers)
        
        
        #pdb.set_trace()
        if n_iterations < n_min:
            converge = False
        if n_iterations> n_min and delta < 1e-6:
            converge = True
            print 'convergence reached\ndelta = {}'.format(delta)
        n_iterations+=1
        
    return Planet
    
    
def initialize_by_mass(*args):
    
    M_P                      = args[0]*Earth_mass
    structural_params        = args[1]
    compositional_params     = args[2]
    num_mantle_layers, num_core_layers, number_h2o_layers = args[3]
    Mf_core = args[4] 

    wt_frac_water, FeMg, SiMg, CaMg, AlMg, mol_frac_Fe_mantle, wt_frac_Si_core, \
     wt_frac_O_core, wt_frac_S_core = compositional_params

    core_rad_frac         = structural_params[6]
    Mantle_potential_temp = structural_params[7]
    water_rad_frac        = structural_params[8]
    water_potential_temp  = structural_params[9]

    #mass of each layer
    M_h20    = M_P*wt_frac_water
    M_core   = M_P*Mf_core
    M_mantle = M_P-(M_h20+M_core)
    
    #total number of layers 
    n_tot = num_core_layers+num_mantle_layers + number_h2o_layers
    num_layers = n_tot
    
    
    # these grids are initialized in this function and are then a passed around through the other routines and edited each iteration
    radius_layers      = np.zeros(num_layers)
    density_layers     = np.zeros(num_layers)
    volume_layers      = np.zeros(num_layers)
    mass_layers        = np.zeros(num_layers)
    cumulative_mass    = np.zeros(num_layers)
    Temperature_layers = np.zeros(num_layers)
    
    # used for compressiofh2on funciton
    gravity_layers  = np.zeros(num_layers)
    Pressure_layers = np.zeros(num_layers)
    alpha           = np.zeros(num_layers)
    cp              = np.zeros(num_layers)
    Vphi            = np.zeros(num_layers)

    Vs = np.zeros(num_layers)
    Vp = np.zeros(num_layers)
    K  =  np.zeros(num_layers)

    #initialize density, volume and radius of core, mantle and water layer
    rhoCor     = 8280.0
    V_core   = M_core / rhoCor
    R_core   = ((0.75 * V_core) / np.pi) ** (1.0 / 3.0)

    # Mantle
    rhoMan = 4100.
    V_man  = (M_mantle / rhoMan)
    R_man  = (((0.75 * V_man) / np.pi) + R_core ** 3) ** (1.0 / 3.0) #radius up to top of mantle

    # Water
    rhoH20       = 1100.
    V_h2o        = M_h20/rhoH20
    R_h2o        = (((0.75 * V_h2o) / np.pi) + R_man ** 3) ** (1.0 / 3.0)
    R_P          = R_h2o #radius of the planet 


    #planet_radius_guess    = radius_planet*Earth_radius
    #water_thickness_guess  = water_rad_frac*planet_radius_guess
    #core_thickness_guess   = core_rad_frac * (planet_radius_guess-water_thickness_guess)
    #mantle_thickness_guess = planet_radius_guess - water_thickness_guess - core_thickness_guess



    #define center of planet
    radius_layers[0]      = 0.0
    density_layers[0]     = 8280.0
    Temperature_layers[0] = 4500.
    volume_layers[0]      = (4.0 / 3) * np.pi * (radius_layers[0] ** 3)
    mass_layers[0]        = volume_layers[0]*density_layers[0]
    cumulative_mass[0]    = mass_layers[0]





    for i in range(1, n_tot):
        
        if i <num_core_layers:
            radius_layers[i]      = ((float(i)/(num_core_layers-1))*R_core)
            density_layers[i]     = rhoCor
            Temperature_layers[i] = 4500.+ (float((3000.-4500.))/float(num_core_layers-1.))\
                                                            *float((i))
            volume_layers[i]      = (4.0 / 3) * np.pi * (radius_layers[i] ** 3 - radius_layers[i-1] ** 3)
            mass_layers[i]        = volume_layers[i]*density_layers[i]
            cumulative_mass[i]    = mass_layers[i]+cumulative_mass[i-1]
            
        elif i < (num_core_layers+num_mantle_layers):
            #print (float(i-num_core_layers)/(num_mantle_layers-1))
            radius_layers[i]      = (R_core+((float(i-num_core_layers)/(num_mantle_layers-1.))*(R_man-R_core)))
            density_layers[i]     = rhoMan
            Temperature_layers[i] = 2700.
            volume_layers[i]      = (4.0 / 3) * np.pi * (radius_layers[i] ** 3 - radius_layers[i-1] ** 3)
            mass_layers[i]        = volume_layers[i]*density_layers[i]
            cumulative_mass[i]    = mass_layers[i]+cumulative_mass[i-1]

        elif number_h2o_layers>0:
            radius_layers[i]      = R_man+\
                                   ((float(i-num_core_layers-num_mantle_layers)/(number_h2o_layers-1.))*(R_h2o-R_man))
            density_layers[i]     = rhoH20
            Temperature_layers[i] = 300.
            volume_layers[i]      = (4.0 / 3) * np.pi * (radius_layers[i] ** 3 - radius_layers[i-1] ** 3)
            mass_layers[i]        = volume_layers[i]*density_layers[i]
            cumulative_mass[i]    = mass_layers[i]+cumulative_mass[i-1]


    for i in range(num_layers):

        if i > num_core_layers+num_mantle_layers:
            Pressure_layers[i] = 1
        else:
            Pressure_layers[i] = (float((5000.-(300.*10000))/float(num_core_layers+num_mantle_layers-1))*float(i)
                                  + 300.*10000)
            #print 'P[{}] = {}'.format(i,Pressure_layers[i])

    keys = ['radius','density','temperature','gravity','pressure',\
            'alpha','cp','Vphi','Vp','Vs','K', 'dmass', 'mass', 'volume']


                          
    return dict(zip(keys,[radius_layers, density_layers,Temperature_layers,gravity_layers, Pressure_layers,
                          alpha, cp,Vphi,Vp,Vs,K, mass_layers, cumulative_mass, volume_layers]))




