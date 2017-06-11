
import numpy as np
import burnman
from scipy import interpolate
import math
from scipy.integrate import odeint
import sys


ToPa = 100000.
ToBar = 1./ToPa
G = 6.67408e-11

def get_rho(Planet,grids,Core_wt_per,structural_params):

    Pressure_layers = Planet.get('pressure')
    Temperature_layers = Planet.get('temperature')
    rho_layers = Planet.get('density')

    num_mantle_layers, num_core_layers, number_h2o_layers = structural_params[4:]

    num_layers = num_mantle_layers+num_core_layers+number_h2o_layers

    P_points_mantle = np.zeros(num_mantle_layers)
    T_points_mantle = np.zeros(num_mantle_layers)

    for i in range(num_core_layers):
        if i <= num_core_layers:
            rho_layers[i] = get_core_rho(Pressure_layers[i],Temperature_layers[i],Core_wt_per)

    for i in range(num_mantle_layers):
        P_points_mantle[i] = Pressure_layers[i+num_core_layers]
        T_points_mantle[i] = Temperature_layers[i+num_core_layers]

    mantle_data = interpolate.griddata((grids['pressure'], grids['temperature']),
                                           grids['density'],(P_points_mantle, T_points_mantle), method='linear')

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

def get_water_rho(Pressure,Temperature):
    return 0

def get_gravity(Planet):
    radii = Planet.get('radius')
    density = Planet.get('density')
    rhofunc = interpolate.UnivariateSpline(radii, density)
    # Create a spline fit of density as a function of radius

    # Numerically integrate Poisson's equation
    poisson = lambda p, x: 4.0 * np.pi * G * rhofunc(x) * x * x

    gravity_layers = np.ravel(odeint(poisson, 0., radii))
    gravity_layers[1:] = gravity_layers[1:]/radii[1:]/radii[1:]
    gravity_layers[0] = 0

    return gravity_layers

def get_pressure(Planet):
    radii = Planet.get('radius')
    density = Planet.get('density')
    gravity = Planet.get('gravity')

    density = density
    # convert radii to depths
    depths = radii[-1] - radii
    # Make a spline fit of density as a function of depth
    rhofunc = interpolate.UnivariateSpline(depths[::-1], density[::-1])
    # Make a spline fit of gravity as a function of depth
    gfunc = interpolate.UnivariateSpline(depths[::-1], gravity[::-1])

    # integrate the hydrostatic equation
    pressure = np.ravel(odeint((lambda p, x: gfunc(x) * rhofunc(x)),1.e9, depths[::-1]))

    pressure = [((i/1.e9)*10000.) for i in pressure]
    return np.asarray(pressure[::-1])

def get_mass(Planet):
    radii = Planet.get('radius')
    density = Planet.get('density')
    rhofunc = interpolate.UnivariateSpline(radii, density)
    mass_in_sphere = lambda p, x: 4.0 * np.pi * rhofunc(x) * x * x

    mass = np.ravel(odeint(mass_in_sphere,0.,radii))
    return mass

def get_temperature(Planet,grids,structural_params):
    radii = Planet.get('radius')
    gravity = Planet.get('gravity')
    temperature = Planet.get('temperature')
    pressure = Planet.get('pressure')

    Mantle_potential_temp,num_mantle_layers, num_core_layers, number_h2o_layers = structural_params[3:]

    radii = radii[num_core_layers:]
    gravity =  gravity[num_core_layers:]

    pressure = pressure[num_core_layers:]

    temperature = temperature[num_core_layers:]
    depths = radii[-1] - radii

    spec_heat= interpolate.griddata((grids['pressure'],grids['temperature']),grids['cp'],
                                    (pressure,temperature),method='linear')
    alpha= interpolate.griddata((grids['pressure'],grids['temperature']),grids['alpha'],
                                (pressure,temperature),method='linear')


    grav_func = interpolate.UnivariateSpline(depths[::-1],gravity[::-1])
    spec_heat_func = interpolate.UnivariateSpline(depths[::-1],spec_heat[::-1])
    alpha_func = interpolate.UnivariateSpline(depths[::-1],alpha[::-1])

    adiabat = lambda p, x:  alpha_func(x)*grav_func(x) / spec_heat_func(x)
    initial_grad = alpha[0]*gravity[0]/spec_heat[0]

    gradient = np.ravel(odeint(adiabat, initial_grad,depths[::-1]))

    mantle_temperatures = [math.exp(i)*Mantle_potential_temp for i in gradient][::-1]

    core_temperatures = Planet['temperature'][:num_core_layers]
    temperatures = np.append(core_temperatures,mantle_temperatures)

    return temperatures


def check_convergence(new_rho,old_rho):

    delta = ([(1.-(old_rho[i]/new_rho[i])) for i in range(len(new_rho))])
    new_rho = [i for i in new_rho]

    for i in delta:
        if i >= 1.e-6:
            print i
            return False,new_rho
        else:
            return True, new_rho

