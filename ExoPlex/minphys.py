
import numpy as np
import burnman
from scipy import interpolate
import math
from scipy.integrate import odeint
import sys
import functions
import pdb


ToPa = 100000.
ToBar = 1./ToPa
G = 6.67408e-11

def get_rho(Planet,grids,Core_wt_per,layers):

    Pressure_layers    = Planet.get('pressure')
    Temperature_layers = Planet.get('temperature')
    rho_layers         = Planet.copy().get('density')

    num_mantle_layers, num_core_layers, number_h2o_layers = layers

    num_layers = num_mantle_layers+num_core_layers+number_h2o_layers

    #P_points_mantle = np.zeros(num_mantle_layers)
    #T_points_mantle = np.zeros(num_mantle_layers)

    for i in range(num_core_layers):
        if i <= num_core_layers:
            rho_layers[i] = get_core_rho(Pressure_layers[i],Temperature_layers[i],Core_wt_per)

    if number_h2o_layers > 0:
        P_points_water = Pressure_layers[(num_mantle_layers+num_core_layers):]
        T_points_water = Temperature_layers[(num_mantle_layers+num_core_layers):]


        water_rho = get_water_rho(P_points_water, T_points_water)
        
        rho_layers[num_core_layers+num_mantle_layers:] = water_rho

    P_points_UM = []
    T_points_UM = []
    P_points_LM = []
    T_points_LM = []

    for i in range(num_mantle_layers):
        #DEBUG
        #print 'P[{:.0f}] = {:.3f}]'.format(i+num_core_layers,Pressure_layers[i+num_core_layers])
        if Pressure_layers[i+num_core_layers] >=1250000:
            P_points_LM.append(Pressure_layers[i+num_core_layers])
            T_points_LM.append(Temperature_layers[i+num_core_layers])
        else:
            P_points_UM.append(Pressure_layers[i+num_core_layers])
            T_points_UM.append(Temperature_layers[i+num_core_layers])

    UM_data = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                           grids[0]['density'],(P_points_UM, T_points_UM), method='linear')

    LM_data = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                           grids[1]['density'],(P_points_LM, T_points_LM), method='linear')

    to_switch_P = []
    to_switch_T = []
    to_switch_index = []
    for i in range(len(UM_data)):
        if np.isnan(UM_data[i]) == True:
            #try lower mantle:
            to_switch_P.append(P_points_UM[i])
            to_switch_T.append(T_points_UM[i])
            to_switch_index.append(i)

    if len(to_switch_P) >0:
        test = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                    grids[1]['density'], (to_switch_P, to_switch_T), method='linear')

        for i in range(len(test)):

            if np.isnan(test[i]) == True:
                print to_switch_P[i], to_switch_T[i]
                print "UM Rho Outside of range!"
                sys.exit()
            else:
                UM_data[to_switch_index[i]] = test[i]

    to_switch_P = []
    to_switch_T = []
    to_switch_index = []

    for i in range(len(LM_data)):
        if np.isnan(LM_data[i]) == True:
            # try upper mantle:
            to_switch_P.append(P_points_LM[i])
            to_switch_T.append(T_points_LM[i])
            to_switch_index.append(i)

    if len(to_switch_P) > 0:
        test = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                    grids[0]['density'], (to_switch_P, to_switch_T), method='linear')

        for i in range(len(test)):

            if np.isnan(test[i]) == True:
                print to_switch_P[i], to_switch_T[i]
                print to_switch_T
                print '\n\n'
                print min(Pressure_layers)
                print "LM Rho Outside of range!"
                sys.exit()
            else:
                LM_data[to_switch_index[i]] = test[i]

    mantle_data = np.append(LM_data,UM_data)

    rho_layers[num_core_layers:num_core_layers+num_mantle_layers] = mantle_data


    return rho_layers

def get_water_rho(Pressure,Temperature):
    phase = []
    P_water = []
    T_water = []
    P_ice = []
    T_ice = []
    density_ice = []
    for i in range(len(Pressure)):
        if (functions.find_water_phase(Pressure[i],Temperature[i]) == 'Water'):
            P_water.append(Pressure[i]*(1.e9/10000.)-((73.5e-5)*2.06e9*(Temperature[i]-373.)))
            T_water.append(Temperature[i])

        else:
            phase.append(functions.find_water_phase(Pressure[i],Temperature[i]))
            P_ice.append(Pressure[i])
            T_ice.append(Temperature[i])

    #cols = np.arange(1, 701)
    #liqDat = np.genfromtxt('../ExoPlex/water_table.dat', usecols=cols)
    #Twat = np.arange(300, 1000)
    #Pwat = np.genfromtxt('../ExoPlex/water_table.dat', usecols=(0,))
    #tab = np.zeros((len(Twat)*len(Pwat),3))
    #for i in range(len(Twat)):

    #for j in range(len(Pwat)):
    #        tab[((i*len(Pwat)))+j] = [Pwat[j],Twat[i],liqDat[j][i]]

    rock=burnman.minerals.other.water()
    density_water = rock.evaluate(['density'],P_water,T_water)[0]
    for i in range(len(phase)):
        if phase[i] == 'Ice_VII':
            class Ice_VII(burnman.Mineral):

                def __init__(self):
                    self.params = {
                        'name': 'ice_VII',
                        'equation_of_state': 'bm2',
                        'V_0':12.49e-6,
                        'K_0': 20.15e9,
                        'Kprime_0':4.,
                        'molar_mass': 0.01801528,
                    }
                    burnman.Mineral.__init__(self)

            rock = Ice_VII()
            density_ice.append(rock.evaluate(['density'], 1.e9 * (P_ice[i] / 10000.), T_ice[i])[0])

        if phase[i] == 'Ice_VI':
            class Ice_VI(burnman.Mineral):
                def __init__(self):
                    self.params = {
                        'name': 'ice_VI',
                        'equation_of_state': 'bm2',
                        'V_0': 14.17e-6,
                        'K_0': 14.01e9,
                        'Kprime_0': 4.,
                        'molar_mass': 0.01801528,
                    }
                    burnman.Mineral.__init__(self)

            rock = Ice_VI()
            density_ice.append(rock.evaluate(['density'], 1.e9 * (P_ice[i] / 10000.), T_ice[i])[0])

        if phase[i] == 'Ice_Ih':
            class Ice_Ih(burnman.Mineral):

                def __init__(self):
                    self.params = {
                        'name': 'Ice_Ih',
                        'equation_of_state': 'bm3',
                        'V_0': 1./(916.72/0.01801528),
                        'K_0': 9.2e9,
                        'Kprime_0': 5.5,
                        'molar_mass': 0.01801528,
                    }
                    burnman.Mineral.__init__(self)

            rock = Ice_Ih()
            print "uh oh"
            density_ice.append(rock.evaluate(['density'], 1.e9 * (P_ice[i] / 10000.), T_ice[i])[0])

    density = np.concatenate((density_ice,density_water),axis=0)
    return density

def get_water_Cp(Pressure, Temperature):
    phase = functions.find_water_phase(Pressure,Temperature)
    Pressure = Pressure/10000.
    if phase == 'Water':
        return 4.184e3

    if phase == 'Ice_VII' or phase=='Ice_VI':
        cp = 3.3 + 22.1 * np.exp(-0.058 * Pressure)  # Asahara 2010
        return 1000.*cp

    if phase == 'Ice_Ih':
        return 4.184e3


def get_water_alpha(Pressure,Temperature):
    phase = functions.find_water_phase(Pressure, Temperature)
    Pressure = Pressure / 10000.
    if phase == 'Water':
        return 214.e-6

    if phase == 'Ice_VII' or phase=='Ice_VI':
        Ks = 23.7
        Ksp = 4.15
        a0 = -3.9e-4
        a1 = 1.5e-6
        at = a0 + a1 * Temperature
        alpha = at * (1 + (Ksp / Ks) * Pressure) ** (-0.9)  # fei 1993
        return alpha

    if phase == 'Ice_Ih':
        return 214.e-6

    return 0

def get_core_rho(Pressure,Temperature,Core_wt_per):
    wt_frac_Si = Core_wt_per.get('Si')
    wt_frac_O = Core_wt_per.get('O')
    wt_frac_S = Core_wt_per.get('S')
    wt_frac_Fe = Core_wt_per.get('Fe')

    mFe = 55.845 #molar weights
    mSi = 28.0867
    mO = 15.9994
    mS = 32.0650

    mol_total = (wt_frac_Fe/mFe)+(wt_frac_O/mO)+(wt_frac_S/mS)+(wt_frac_Si/mSi)
    mol_frac_Fe = (wt_frac_Fe/mFe) / mol_total

    mol_frac_Si = (wt_frac_Si/mSi) / mol_total
    mol_frac_S = (wt_frac_S/mS) / mol_total
    mol_frac_O = (wt_frac_O/mO) / mol_total

    molar_weight_core = (mol_frac_Fe*mFe) + (mol_frac_Si * mSi) + (mol_frac_O*mO) + (mol_frac_S*mS)

    class iron(burnman.Mineral):

        def __init__(self):
            self.params = {
                'name': 'iron',
                'equation_of_state': 'bm4',
                'V_0': 7.95626e-6,
                'K_0': 109.7e9,
                'Kprime_0': 4.66,
                'Kprime_prime_0': -0.043e-9,
                'molar_mass': molar_weight_core/1000.,
            }
            burnman.Mineral.__init__(self)

    rock = iron()
    density = rock.evaluate(['density'], 1.e9*(Pressure/10000.), Temperature)[0]
    return density

def get_core_speeds(Pressure,Temperature,Core_wt_per):
    wt_frac_Si = Core_wt_per.get('Si')
    wt_frac_O = Core_wt_per.get('O')
    wt_frac_S = Core_wt_per.get('S')
    wt_frac_Fe = Core_wt_per.get('Fe')

    mFe = 55.845 #molar weights
    mSi = 28.0867
    mO = 15.9994
    mS = 32.0650

    mol_total = (wt_frac_Fe/mFe)+(wt_frac_O/mO)+(wt_frac_S/mS)+(wt_frac_Si/mSi)
    mol_frac_Fe = (wt_frac_Fe/mFe) / mol_total

    mol_frac_Si = (wt_frac_Si/mSi) / mol_total
    mol_frac_S = (wt_frac_S/mS) / mol_total
    mol_frac_O = (wt_frac_O/mO) / mol_total

    molar_weight_core = (mol_frac_Fe*mFe) + (mol_frac_Si * mSi) + (mol_frac_O*mO) + (mol_frac_S*mS)

    class iron(burnman.Mineral):

        def __init__(self):
            self.params = {
                'name': 'iron',
                'equation_of_state': 'bm4',
                'V_0': 7.95626e-6,
                'K_0': 109.7e9,
                'Kprime_0': 4.66,
                'Kprime_prime_0': -0.043e-9,
                'molar_mass': molar_weight_core/1000.,
            }
            burnman.Mineral.__init__(self)

    Pressure = [i*((1.e9)/10000.) for i in Pressure]
    Temperature = [i for i in Temperature]

    rock = iron()
    speeds = rock.evaluate(['v_phi', 'v_p', 'v_s'], Pressure, Temperature)

    return speeds

def calc_gravity(Planet, layers):
    mass    = Planet.get('mass') #cumulative mass
    rad     = Planet.get('radius') 
    grav    = Planet.get('gravity')

    num_mantle_layers, num_core_layers, number_h2o_layers = layers
    
    n_tot = num_mantle_layers + num_core_layers + number_h2o_layers
    
    # gravity at each zone
    for i in range(1, n_tot):
        grav[i] = (G * mass[i]) / (rad[i] ** 2)

    return grav
    

def get_gravity(Planet,layers):
    radii   = Planet.get('radius')
    density = Planet.get('density')


    num_mantle_layers, num_core_layers, number_h2o_layers = layers

    
    
    radii_core =  radii[:num_core_layers]
    
    density_core = density[:num_core_layers]

    radii_mantle = radii[num_core_layers:(num_core_layers+num_mantle_layers)]
    density_mantle = density[num_core_layers:(num_core_layers+num_mantle_layers)]

    radii_water = radii [(num_core_layers+num_mantle_layers):]
    density_water = density[(num_core_layers+num_mantle_layers):]

    rhofunc_core = interpolate.UnivariateSpline(radii_core, density_core)
    rhofunc_mantle = interpolate.UnivariateSpline(radii_mantle, density_mantle)


    # Create a spline fit of density as a function of radius

    # Numerically integrate Poisson's equation
    poisson_core   = lambda p, x: 4.0 * np.pi * G * rhofunc_core(x) * x * x
    poisson_mantle = lambda p, x: 4.0 * np.pi * G * rhofunc_mantle(x) * x * x

    gravity_layers_core = np.ravel(odeint(poisson_core, 0., radii_core))
    gravity_layers_mantle = np.ravel(odeint(poisson_mantle,gravity_layers_core[-1],radii_mantle))

    if number_h2o_layers > 0:
        rhofunc_water = interpolate.UnivariateSpline(radii_water, density_water)
        poisson_water = lambda p, x: 4.0 * np.pi * G * rhofunc_water(x) * x * x
        gravity_layers_water = np.ravel(odeint(poisson_water,gravity_layers_mantle[-1],radii_water))
        gravity_layers = np.concatenate((gravity_layers_core,gravity_layers_mantle,gravity_layers_water),axis=0)
        #print gravity_layers_water

    else:
        gravity_layers = np.concatenate((gravity_layers_core, gravity_layers_mantle), axis=0)


    gravity_layers[1:] = gravity_layers[1:]/radii[1:]/radii[1:]
    gravity_layers[0] = 0

    return gravity_layers

def calc_pressure(Planet, layers):
    #pressure comes in as bar. 
    #calculate pressure in Pa and convert to bar
    
    radius   = Planet.get('radius')
    density  = Planet.get('density')
    gravity  = Planet.get('gravity')
    pressure = Planet.get('pressure')

    num_mantle_layers, num_core_layers, number_h2o_layers = layers

    n_tot = num_mantle_layers + num_core_layers + number_h2o_layers


    P_pa  = np.zeros(n_tot)
    P_bar = np.zeros(n_tot)
    #keep surface pressure fixed, change to Pa
    P_pa[n_tot-1] = pressure[-1]*ToPa

    for i in range(n_tot - 2, -1, -1):
        gmid   = 0.5*(gravity[i+1]+gravity[i])
        rhomid = 0.5*(density[i+1]+density[i])
        P_pa[i] = P_pa[i + 1] + gmid*rhomid * (radius[i + 1] - radius[i]) 
        
    P_bar = P_pa*ToBar
    
    return P_bar
def get_pressure(Planet,layers):
    radii   = Planet.get('radius')
    density = Planet.get('density')
    gravity = Planet.get('gravity')

    num_mantle_layers, num_core_layers, number_h2o_layers = layers


    # convert radii to depths
    depths       = radii[-1] - radii
    depths_core  = depths[:num_core_layers]
    gravity_core = gravity[:num_core_layers]
    density_core = density[:num_core_layers]

    depths_mant  = depths[num_core_layers:(num_core_layers+num_mantle_layers)]
    gravity_mant = gravity[num_core_layers:(num_core_layers+num_mantle_layers)]
    density_mant = density[num_core_layers:(num_core_layers+num_mantle_layers)]



    # Make a spline fit of density as a function of depth
    rhofunc_mant = interpolate.UnivariateSpline(depths_mant[::-1], density_mant[::-1])
    # Make a spline fit of gravity as a function of depth
    gfunc_mant   = interpolate.UnivariateSpline(depths_mant[::-1], gravity_mant[::-1])

    rhofunc_core = interpolate.UnivariateSpline(depths_core[::-1], density_core[::-1])
    gfunc_core   = interpolate.UnivariateSpline(depths_core[::-1], gravity_core[::-1])

    if number_h2o_layers > 0:
        depths_water  = depths[(num_core_layers + num_mantle_layers):]
        gravity_water = gravity[(num_core_layers + num_mantle_layers):]
        density_water = density[(num_core_layers + num_mantle_layers):]

        rhofunc_water = interpolate.UnivariateSpline(depths_water[::-1], density_water[::-1])
        gfunc_water   = interpolate.UnivariateSpline(depths_water[::-1], gravity_water[::-1])

        #integrate from 1 bar
        pressure_water = np.ravel(odeint((lambda p, x: gfunc_water(x) * rhofunc_water(x)),(1./10000.)*1.e9, depths_water[::-1]))
        WMB_pres       = pressure_water[-1] #does not account for very small water layers and ends up breaking code
        WMB_pres = 5.e8


    else:
        WMB_pres = 5.e8



    # integrate the hydrostatic equation
    pressure_mant = np.ravel(odeint((lambda p, x: gfunc_mant(x) * rhofunc_mant(x)),WMB_pres, depths_mant[::-1]))
    CMB_pres      = pressure_mant[-1]
    pressure_core = np.ravel(odeint((lambda p, x: gfunc_core(x) * rhofunc_core(x)),CMB_pres, depths_core[::-1]))

    if number_h2o_layers > 0:
        pressure = np.concatenate((pressure_water,pressure_mant,pressure_core),axis=0)
        pressure = [((i/1.e9)*10000.) for i in pressure]
        return np.asarray(pressure[::-1])
    else:
        pressure = np.concatenate((pressure_mant,pressure_core),axis=0)
        pressure = [((i/1.e9)*10000.) for i in pressure]
        return np.asarray(pressure[::-1])


def get_mass(Planet,layers):
    radii   = Planet.get('radius')
    density = Planet.get('density')
    num_mantle_layers, num_core_layers, number_h2o_layers = layers

    if len(radii)<=num_core_layers:

        rhofunc = interpolate.UnivariateSpline(radii, density)
        mass_in_sphere = lambda p, x: 4.0 * np.pi * rhofunc(x) * x * x

        mass = np.ravel(odeint(mass_in_sphere,0.,radii))
        return mass

    elif len(radii) > (num_core_layers+num_mantle_layers):
        radii_core = radii[:num_core_layers]
        density_core = density[:num_core_layers]
        rhofunc_core = interpolate.UnivariateSpline(radii_core, density_core)
        mass_in_core = lambda p, x: 4.0 * np.pi * rhofunc_core(x) * x * x

        mass_core = np.ravel(odeint(mass_in_core, 0., radii_core))

        radii_mantle = radii[num_core_layers:num_core_layers + num_mantle_layers]
        density_mantle = density[num_core_layers:num_core_layers + num_mantle_layers]
        rhofunc_mantle = interpolate.UnivariateSpline(radii_mantle, density_mantle)
        mass_in_mantle = lambda p, x: 4.0 * np.pi * rhofunc_mantle(x) * x * x

        mass_mantle = np.ravel(odeint(mass_in_mantle,mass_core[-1],radii_mantle))

        radii_water = radii[num_core_layers + num_mantle_layers:]
        density_water = density[num_core_layers + num_mantle_layers:]
        rhofunc_water = interpolate.UnivariateSpline(radii_water, density_water)

        mass_in_water = lambda p, x: 4.0 * np.pi * rhofunc_water(x) * x * x

        mass_water= np.ravel(odeint(mass_in_water, mass_mantle[-1], radii_water))

        return np.concatenate((mass_core,mass_mantle,mass_water),axis=0)
    else:
        radii_core = radii[:num_core_layers]
        density_core = density[:num_core_layers]
        rhofunc_core = interpolate.UnivariateSpline(radii_core, density_core)
        mass_in_core = lambda p, x: 4.0 * np.pi * rhofunc_core(x) * x * x

        mass_core = np.ravel(odeint(mass_in_core, 0., radii_core))

        radii_mantle = radii[num_core_layers:(num_core_layers + num_mantle_layers)]
        density_mantle = density[num_core_layers:(num_core_layers + num_mantle_layers)]
        rhofunc_mantle = interpolate.UnivariateSpline(radii_mantle, density_mantle)
        mass_in_mantle = lambda p, x: 4.0 * np.pi * rhofunc_mantle(x) * x * x

        mass_mantle = np.ravel(odeint(mass_in_mantle,mass_core[-1],radii_mantle))

        return np.concatenate((mass_core,mass_mantle),axis=0)

def get_temperature(Planet,grids,structural_parameters,layers):
    radii       = Planet.get('radius')
    gravity     = Planet.get('gravity')
    temperature = Planet.get('temperature')
    pressure    = Planet.get('pressure')
    num_mantle_layers, num_core_layers, number_h2o_layers = layers

    
    Mantle_potential_temp = structural_parameters[7]
    Water_potential_temp  = structural_parameters[9]

    #temperature gradient is an adiabat in the mantle and water layers, below we seperate the mantle data
    radii       = radii[num_core_layers:]
    gravity     = gravity[num_core_layers:]
    pressure    = pressure[num_core_layers:]
    temperature = temperature[num_core_layers:]
    depths      = radii[-1] - radii

    P_points_UM = []
    T_points_UM = []
    P_points_LM = []
    T_points_LM = []

    P_points_water = pressure[num_mantle_layers:]
    T_points_water = temperature[num_mantle_layers:]

    spec_heat_water = []
    alpha_water     = []

    #find Cp and alpha of water layer
    for i in range(len(P_points_water)):
        spec_heat_water.append(get_water_Cp(P_points_water[i],T_points_water[i]))
        alpha_water.append(get_water_alpha(P_points_water[i],T_points_water[i]))
    
    #find upper and lower mantle points
    for i in range(num_mantle_layers):
        if pressure[i] >=1250000:
            P_points_LM.append(pressure[i])
            T_points_LM.append(temperature[i])
        else:
            P_points_UM.append(pressure[i])
            T_points_UM.append(temperature[i])


    depths_mantle  = depths[:num_mantle_layers]
    gravity_mantle = gravity[:num_mantle_layers]

    depths_water  = depths[num_mantle_layers:]
    gravity_water = gravity[num_mantle_layers:]

    UM_cp_data = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                           grids[0]['cp'],(P_points_UM, T_points_UM), method='linear')

    LM_cp_data = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                           grids[1]['cp'],(P_points_LM, T_points_LM), method='linear')

    to_switch_P = []
    to_switch_T = []
    to_switch_index = []
    for i in range(len(UM_cp_data)):
        if np.isnan(UM_cp_data[i]) == True:
            # try lower mantle:
            to_switch_P.append(P_points_UM[i])
            to_switch_T.append(T_points_UM[i])
            to_switch_index.append(i)


    if len(to_switch_P) > 0:
        test = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                    grids[1]['cp'], (to_switch_P, to_switch_T), method='linear')


        for i in range(len(test)):

            if np.isnan(test[i]) == True:
                print to_switch_P[i] , to_switch_T[i]
                print "UM Cp Outside of range!"
                sys.exit()
            else:
                UM_cp_data[to_switch_index[i]] = test[i]

    to_switch_P = []
    to_switch_T = []
    to_switch_index = []

    for i in range(len(LM_cp_data)):
        if np.isnan(LM_cp_data[i]) == True:
            # try lower mantle:
            to_switch_P.append(P_points_LM[i])
            to_switch_T.append(T_points_LM[i])
            to_switch_index.append(i)

    if len(to_switch_P) > 0:
        test = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                    grids[0]['cp'], (to_switch_P, to_switch_T), method='linear')

        for i in range(len(test)):
            if np.isnan(test[i]) == True:
                print to_switch_P[i], to_switch_T[i]
                print "LM Cp Outside of range!"
                sys.exit()
            else:
                LM_cp_data[to_switch_index[i]] = test[i]

    spec_heat_mantle = np.append(LM_cp_data,UM_cp_data)

    UM_alpha_data = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                           grids[0]['alpha'],(P_points_UM, T_points_UM), method='linear')

    LM_alpha_data = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                           grids[1]['alpha'],(P_points_LM, T_points_LM), method='linear')

    to_switch_P = []
    to_switch_T = []
    to_switch_index = []
    for i in range(len(UM_alpha_data)):
        if np.isnan(UM_alpha_data[i]) == True:
            # try lower mantle:
            to_switch_P.append(P_points_UM[i])
            to_switch_T.append(T_points_UM[i])
            to_switch_index.append(i)

    if len(to_switch_P) > 0:
        test = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                    grids[1]['alpha'], (to_switch_P, to_switch_T), method='linear')

        for i in range(len(test)):

            if np.isnan(test[i]) == True:
                print to_switch_P[i] / 1e5, to_switch_T[i]
                print "UM Alpha Outside of range!"
                sys.exit()
            else:
                UM_alpha_data[to_switch_index[i]] = test[i]

    to_switch_P = []
    to_switch_T = []
    to_switch_index = []
    for i in range(len(LM_alpha_data)):
        if np.isnan(LM_alpha_data[i]) == True:
            # try lower mantle:
            to_switch_P.append(P_points_LM[i])
            to_switch_T.append(T_points_LM[i])
            to_switch_index.append(i)

    if len(to_switch_P) > 0:
        test = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                    grids[0]['alpha'], (to_switch_P, to_switch_T), method='linear')

        for i in range(len(test)):

            if np.isnan(test[i]) == True:
                print to_switch_P[i] / 1e5, to_switch_T[i]
                print "LM Alpha Outside of range!"
                sys.exit()
            else:
                LM_alpha_data[to_switch_index[i]] = test[i]

    alpha_mantle = np.append(LM_alpha_data,UM_alpha_data)

    grav_func      = interpolate.UnivariateSpline(depths_mantle[::-1],gravity_mantle[::-1])
    spec_heat_func = interpolate.UnivariateSpline(depths_mantle[::-1],spec_heat_mantle[::-1])
    alpha_func     = interpolate.UnivariateSpline(depths_mantle[::-1],alpha_mantle[::-1])

    adiabat_mantle = lambda p, x:  alpha_func(x)*grav_func(x) / spec_heat_func(x)

    gradient_mantle = np.ravel(odeint(adiabat_mantle, 0.,depths_mantle[::-1]))


    """
    for i in range(len(depths_mantle)):
        if np.isnan(spec_heat_mantle[i]) == True:
            print "There's a nan in Cp"
            print i, pressure[i]/(1e5), temperature[i],depths[i]/depths[0]
            print "pressure range mantle", structural_parameters[0]
            print "temperature range mantle", structural_parameters[1]
#            print "LM", LM_cp_data
#            print "UM", UM_cp_data
            print len(LM_cp_data)+len(UM_cp_data)
            for j in range(len(LM_cp_data)):
                print P_points_LM[j]/1e5,LM_cp_data[j]
            import matplotlib.pyplot as plt
            plt.plot(pressure,temperature)
            plt.show()
            sys.exit()
    """
    mantle_temperatures = [math.exp(k)*Mantle_potential_temp for k in gradient_mantle][::-1]

    core_temperatures = Planet['temperature'][:num_core_layers]

    if number_h2o_layers > 0:
        grav_func = interpolate.UnivariateSpline(depths_water[::-1], gravity_water[::-1])
        spec_heat_func = interpolate.UnivariateSpline(depths_water[::-1], spec_heat_water[::-1])
        alpha_func = interpolate.UnivariateSpline(depths_water[::-1], alpha_water[::-1])

        adiabat_water = lambda p, x: alpha_func(x) * grav_func(x) / spec_heat_func(x)

        gradient_water = np.ravel(odeint(adiabat_water, 0., depths_water[::-1]))
        water_temperatures = [math.exp(i)*Water_potential_temp for i in gradient_water][::-1]

        return np.concatenate((core_temperatures,mantle_temperatures,water_temperatures),axis=0)

    else:
        return np.concatenate((core_temperatures, mantle_temperatures),axis=0)

def check_convergence(new_rho,old_rho):

    delta = ([(1.-(old_rho[i]/new_rho[i])) for i in range(len(new_rho))])
    new_rho = [i for i in new_rho]

    for i in delta:
        if i >= 1.e-6:
            return False,new_rho
        else:
            return True, new_rho

