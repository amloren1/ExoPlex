
import numpy as np
import burnman
from scipy import interpolate
import math
from scipy.integrate import odeint
import sys


ToPa = 100000.
ToBar = 1./ToPa
G = 6.67408e-11

def get_rho(Planet,grids,Core_wt_per,layers):

    Pressure_layers = Planet.get('pressure')
    Temperature_layers = Planet.get('temperature')
    rho_layers = Planet.get('density')

    num_mantle_layers, num_core_layers, number_h2o_layers = layers

    num_layers = num_mantle_layers+num_core_layers+number_h2o_layers

    #P_points_mantle = np.zeros(num_mantle_layers)
    #T_points_mantle = np.zeros(num_mantle_layers)

    for i in range(num_core_layers):
        if i <= num_core_layers:
            rho_layers[i] = get_core_rho(Pressure_layers[i],Temperature_layers[i],Core_wt_per)

    P_points_UM = []
    T_points_UM = []
    P_points_LM = []
    T_points_LM = []

    for i in range(num_mantle_layers):
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


    mantle_data = np.append(LM_data,UM_data)

    rho_layers[num_core_layers:num_core_layers+num_mantle_layers] = mantle_data
    #else:
    #    rho_layers[i] = get_water_rho(Pressure_layers[i],Temperature_layers[i])
    return rho_layers

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

def get_water_rho(Pressure,Temperature):
    return 0

def get_gravity(Planet,layers):
    radii = Planet.get('radius')
    density = Planet.get('density')
    num_mantle_layers, num_core_layers, number_h2o_layers = layers

    radii_core =  radii[:num_core_layers]
    density_core = density[:num_core_layers]

    radii_mantle = radii[num_core_layers:]
    density_mantle = density[num_core_layers:]

    rhofunc_core = interpolate.UnivariateSpline(radii_core, density_core)
    rhofunc_mantle = interpolate.UnivariateSpline(radii_mantle, density_mantle)
    # Create a spline fit of density as a function of radius

    # Numerically integrate Poisson's equation
    poisson_core = lambda p, x: 4.0 * np.pi * G * rhofunc_core(x) * x * x
    poisson_mantle = lambda p, x: 4.0 * np.pi * G * rhofunc_mantle(x) * x * x

    gravity_layers_core = np.ravel(odeint(poisson_core, 0., radii_core))
    gravity_layers_mantle = np.ravel(odeint(poisson_mantle,gravity_layers_core[-1],radii_mantle))

    gravity_layers = np.concatenate((gravity_layers_core,gravity_layers_mantle),axis=0)
    gravity_layers[1:] = gravity_layers[1:]/radii[1:]/radii[1:]
    gravity_layers[0] = 0

    return gravity_layers

def get_pressure(Planet,layers):
    radii = Planet.get('radius')
    density = Planet.get('density')
    gravity = Planet.get('gravity')
    num_mantle_layers, num_core_layers, number_h2o_layers = layers


    # convert radii to depths
    depths = radii[-1] - radii
    depths_core = depths[:num_core_layers]
    gravity_core = gravity[:num_core_layers]
    density_core = density[:num_core_layers]

    depths_mant = depths[num_core_layers:]
    gravity_mant = gravity[num_core_layers:]
    density_mant = density[num_core_layers:]

    # Make a spline fit of density as a function of depth
    rhofunc_mant = interpolate.UnivariateSpline(depths_mant[::-1], density_mant[::-1])
    # Make a spline fit of gravity as a function of depth
    gfunc_mant = interpolate.UnivariateSpline(depths_mant[::-1], gravity_mant[::-1])

    rhofunc_core = interpolate.UnivariateSpline(depths_core[::-1], density_core[::-1])
    gfunc_core = interpolate.UnivariateSpline(depths_core[::-1], gravity_core[::-1])

    # integrate the hydrostatic equation
    pressure_mant = np.ravel(odeint((lambda p, x: gfunc_mant(x) * rhofunc_mant(x)),5.e8, depths_mant[::-1]))
    CMB_pres = pressure_mant[-1]
    pressure_core = np.ravel(odeint((lambda p, x: gfunc_core(x) * rhofunc_core(x)),CMB_pres, depths_core[::-1]))

    pressure = np.concatenate((pressure_mant,pressure_core),axis=0)
    pressure = [((i/1.e9)*10000.) for i in pressure]
    return np.asarray(pressure[::-1])

def get_mass(Planet):
    radii = Planet.get('radius')
    density = Planet.get('density')
    rhofunc = interpolate.UnivariateSpline(radii, density)
    mass_in_sphere = lambda p, x: 4.0 * np.pi * rhofunc(x) * x * x

    mass = np.ravel(odeint(mass_in_sphere,0.,radii))
    return mass

def get_temperature(Planet,grids,structural_parameters,layers):
    radii = Planet.get('radius')
    gravity = Planet.get('gravity')
    temperature = Planet.get('temperature')
    pressure = Planet.get('pressure')

    num_mantle_layers, num_core_layers, number_h2o_layers = layers

    Mantle_potential_temp = structural_parameters[-1]

    radii = radii[num_core_layers:]
    gravity =  gravity[num_core_layers:]

    pressure = pressure[num_core_layers:]

    temperature = temperature[num_core_layers:]
    depths = radii[-1] - radii

    P_points_UM = []
    T_points_UM = []
    P_points_LM = []
    T_points_LM = []

    for i in range(num_mantle_layers):
        if pressure[i] >=1250000:
            P_points_LM.append(pressure[i])
            T_points_LM.append(temperature[i])
        else:
            P_points_UM.append(pressure[i])
            T_points_UM.append(temperature[i])


    UM_cp_data = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                           grids[0]['cp'],(P_points_UM, T_points_UM), method='linear')

    LM_cp_data = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                           grids[1]['cp'],(P_points_LM, T_points_LM), method='linear')

    spec_heat = np.append(LM_cp_data,UM_cp_data)


    UM_alpha_data = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                           grids[0]['alpha'],(P_points_UM, T_points_UM), method='linear')

    LM_alpha_data = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                           grids[1]['alpha'],(P_points_LM, T_points_LM), method='linear')


    alpha = np.append(LM_alpha_data,UM_alpha_data)

    grav_func = interpolate.UnivariateSpline(depths[::-1],gravity[::-1])
    spec_heat_func = interpolate.UnivariateSpline(depths[::-1],spec_heat[::-1])
    alpha_func = interpolate.UnivariateSpline(depths[::-1],alpha[::-1])

    adiabat = lambda p, x:  alpha_func(x)*grav_func(x) / spec_heat_func(x)

    gradient = np.ravel(odeint(adiabat, 0.,depths[::-1]))

    for i in range(len(depths)):
        if np.isnan(spec_heat[i]) == True:
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

    mantle_temperatures = [math.exp(i)*Mantle_potential_temp for i in gradient][::-1]

    core_temperatures = Planet['temperature'][:num_core_layers]
    temperatures = np.append(core_temperatures,mantle_temperatures)

    return temperatures


def check_convergence(new_rho,old_rho):

    delta = ([(1.-(old_rho[i]/new_rho[i])) for i in range(len(new_rho))])
    new_rho = [i for i in new_rho]

    for i in delta:
        if i >= 1.e-6:
            return False,new_rho
        else:
            return True, new_rho

