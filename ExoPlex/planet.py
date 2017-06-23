import sys
import numpy as np
Earth_radius = 6.371e6
import minphys

def initialize_by_radius(*args):
    radius_planet = args[0]
    structural_params = args[1]
    compositional_params = args[2]
    num_mantle_layers, num_core_layers, number_h2o_layers = args[3]

    wt_frac_water, FeMg, SiMg, CaMg, AlMg, mol_frac_Fe_mantle, wt_frac_Si_core, \
     wt_frac_O_core, wt_frac_S_core = compositional_params

    core_rad_frac,Mantle_potential_temp = structural_params[-2:]


    #def Setup(Mp, masMan, masCor, masH2O, nr, tCrstMan, dtManCor):
    # global Th,Tss
    #    global rad, rho, rho0, vol, dmass, mass, delta,\
    #        T, nic, noc, grav, P, rhon, alpha, cp,phase

    # use inputs to find actual mass of mantle and core

    # if there is a water layer, the imput temperature is lowered because that temperature is for the crustal layer
    # also 50 shells are used for the water layer hence the nh20 vaiable

    if wt_frac_water > 0:

        Surface_temp = 300

    else:
        Surface_temp = Mantle_potential_temp


    if wt_frac_water == 0. and number_h2o_layers > 0:
       print "You have layers of water but no water!"
       number_h2o_layers = 0

    num_layers = num_core_layers+num_mantle_layers + number_h2o_layers  # add 50 shells if there is an h2O layer
    # arrays to be used

    # these grids are initialized in this function and are then a passed around through the other routines and edited each iteration
    radius_layers = np.zeros(num_layers)
    density_layers = np.zeros(num_layers)
    volume_layers = np.zeros(num_layers)
    mass_layers = np.zeros(num_layers)
    cumulative_mass = np.zeros(num_layers)
    Temperature_layers = np.zeros(num_layers)
    ########################################################################################################## Move to compression function and set as outputs
    # used for compressiofh2on funciton
    gravity_layers = np.zeros(num_layers)
    Pressure_layers = np.zeros(num_layers)
    alpha = np.zeros(num_layers)
    cp = np.zeros(num_layers)
    Vphi = np.zeros(num_layers )

    Vs = np.zeros(num_layers )
    Vp = np.zeros(num_layers)
    K =  np.zeros(num_layers)

    # 15 mineral phases + 2ice + liquid water #phasechange
    planet_radius_guess = radius_planet*Earth_radius
    core_thickness_guess = core_rad_frac * planet_radius_guess
    mantle_thickness_guess = planet_radius_guess - (wt_frac_water*planet_radius_guess) - core_thickness_guess

    water_thickness_guess = planet_radius_guess - mantle_thickness_guess - core_thickness_guess


    for i in range(num_layers):
        if i <num_core_layers:
            radius_layers[i]=((float(i)/num_core_layers)*core_thickness_guess)
            density_layers[i]=8280.0
            volume_layers[i]=(4.0 / 3.) * np.pi * (pow(radius_layers[i],3.)- pow(radius_layers[i - 1],3.))
            mass_layers[i]=density_layers[i] * volume_layers[i]
            cumulative_mass[i]=cumulative_mass[i - 1] + mass_layers[i]
            Temperature_layers[i] = 4500.+ (float((3000.-4500.))/float(num_core_layers))\
                                                            *float((i))

        elif i <= (num_core_layers+num_mantle_layers):
            radius_layers[i]=(core_thickness_guess+((float(i-num_core_layers)/num_mantle_layers)*mantle_thickness_guess))
            density_layers[i]=4100.
            volume_layers[i]= (4.0 / 3.) * np.pi * (pow(radius_layers[i],3.)- pow(radius_layers[i - 1],3.))
            mass_layers[i]= density_layers[i] * volume_layers[i]
            cumulative_mass[i]= cumulative_mass[i - 1] + mass_layers[i]

            #Temperature_layers[i] = 4000. + (float((Mantle_potential_temp-4000.))/float(num_mantle_layers))\
            #                                                *float((i-num_core_layers))
            Temperature_layers[i] = 2500.

        else:
            radius_layers[i]=core_thickness_guess+mantle_thickness_guess+\
                             ((float(i-num_core_layers-num_mantle_layers)/number_h2o_layers)*water_thickness_guess)
            density_layers[i]=1000.
            volume_layers[i]=(4.0 / 3.) * np.pi * (pow(radius_layers[i],3.)- pow(radius_layers[i - 1],3.))
            mass_layers[i]=density_layers[i] * volume_layers[i]
            cumulative_mass[i]=cumulative_mass[i - 1] + mass_layers[i]
            Temperature_layers[i] = 300.

    for i in range(num_layers):
        if i > num_core_layers+num_mantle_layers:
            Pressure_layers[i] = 1
        else:
            Pressure_layers[i] = (float((5000.-(300.*10000))/float(num_core_layers+num_mantle_layers))*float(i)
                                  + 300.*10000)


    Pressure_layers[num_core_layers+num_mantle_layers-1] = 10000


    #initial temperature guess of 0.5 K per km
    keys = ['radius','density','temperature','gravity','pressure',\
            'alpha','cp','Vphi''Vp','Vs','K']

    return dict(zip(keys,[radius_layers, density_layers,Temperature_layers,gravity_layers, Pressure_layers,
                          alpha, cp,Vphi,Vp,Vs,K]))



def compress(*args):
    Planet = args[0]
    grids = args[1]
    Core_wt_per = args[2]
    structural_params= args[3]
    layers= args[4]
    n_iterations = 1
    max_iterations = 100


    old_rho = [0  for i in range(len(Planet['density']))]
    converge = False
    print
    while n_iterations <= max_iterations and converge == False:
        print "iteration #",n_iterations
        Planet['density'] = minphys.get_rho(Planet,grids,Core_wt_per,layers)

        for i in range(len(Planet['density'])):
            if np.isnan(Planet['density'][i]) == True:
                print "Density has a nan"
                print i, Planet['pressure'][i],Planet['temperature'][i]
                print "pressure range mantle",structural_params[0]
                print "temperature range mantle",structural_params[1]
                sys.exit()

        Planet['gravity'] = minphys.get_gravity(Planet)

        Planet['pressure'] = minphys.get_pressure(Planet)

        Planet['temperature'] = minphys.get_temperature(Planet,grids,structural_params,layers)
        converge,old_rho = minphys.check_convergence(Planet['density'],old_rho)
        n_iterations+=1

    return Planet
