
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
    temperature_grid, pressure_grid, density_grid, speed_grid, alpha_grid, cp_grid, phases_grid = grids

    for i in range(num_layers+1):
        if i <= num_core_layers:
            rho_layers[i] = get_core_rho(Pressure_layers[i],Temperature_layers[i],Core_wt_per)

        elif i<=num_core_layers+num_mantle_layers:
            rho_layers[i] = get_mantle_rho(Pressure_layers[i],Temperature_layers[i],pressure_grid,\
                                           temperature_grid,density_grid)
        else:
            rho_layers[i] = get_water_rho(Pressure_layers[i],Temperature_layers[i])

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

def get_mantle_rho(Pressure,Temperature,pressure_grid,temperature_grid,density_grid):
    density = interpolate.griddata((pressure_grid,temperature_grid),density_grid,(Pressure,Temperature),method='linear')
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

def get_phases(Planet,grids):

    return 0

def get_temperature(Planet,grids,structural_params):
    radii = Planet.get('radius')
    gravity = Planet.get('gravity')
    temperature = Planet.get('temperature')
    pressure = Planet.get('pressure')

    temperature_grid, pressure_grid, density_grid, speed_grid, alpha_grid, cp_grid, phases_grid = grids
    num_mantle_layers, num_core_layers, number_h2o_layers = structural_params[4:]

    radii = radii[num_core_layers:]
    gravity =  gravity[num_core_layers:]

    spec_heat = np.zeros(len(radii))
    alpha =  np.zeros(len(radii))
    pressure = pressure[num_core_layers:]

    temperature = temperature[num_core_layers:]
    depths = radii[-1] - radii
    for i in range(len(pressure)):
        spec_heat[i] = interpolate.griddata((pressure_grid,temperature_grid),cp_grid,(pressure[i],temperature[i]),method='linear')
        alpha[i] = interpolate.griddata((pressure_grid,temperature_grid),alpha_grid,(pressure[i],temperature[i]),method='linear')


    grav_func = interpolate.UnivariateSpline(depths[::-1],gravity[::-1])
    spec_heat_func = interpolate.UnivariateSpline(depths[::-1],spec_heat[::-1])
    alpha_func = interpolate.UnivariateSpline(depths[::-1],alpha[::-1])

    adiabat = lambda p, x:  alpha_func(x)*grav_func(x) / spec_heat_func(x)
    initial_grad = alpha[0]*gravity[0]/spec_heat[0]

    gradient = np.ravel(odeint(adiabat, initial_grad,depths[::-1]))

    mantle_temperatures = [math.exp(i)*1650. for i in gradient][::-1]

    core_temperatures = Planet.get('temperature')[:num_core_layers]
    temperatures = np.append(core_temperatures,mantle_temperatures)

    return temperatures


def check_convergence(new_rho,old_rho):

    delta = sum([(1.-(old_rho[i]/new_rho[i])) for i in range(len(new_rho))])
    new_rho = [i for i in new_rho]
    print "delta",delta

    return delta,new_rho