
import numpy as np
import burnman
from scipy import interpolate

ToPa = 100000.
def get_rho(Planet,mass_planet,grids,Core_wt_per,structural_params):

    radius_layers, rho_layers, rho0_layers, volume_layers, \
     mass_layers, cumulative_mass, Temperature_layers, \
     gravity, Pressure_layers, rhon, alpha, cp, phase=Planet

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

    for i in range(len(rho_layers)):
        print Pressure_layers[i],rho_layers[i]
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
                'molar_mass': molar_weight_core/100.,
            }
            burnman.Mineral.__init__(self)

    rock = iron()
    density = rock.evaluate(['density'], Pressure*ToPa, Temperature)[0]
    return density

def get_mantle_rho(Pressure,Temperature,pressure_grid,temperature_grid,density_grid):
    density = interpolate.griddata((pressure_grid,temperature_grid),density_grid,(Pressure,Temperature),method='linear')
    return density


def get_water_rho(Pressure,Temperature):
    return 0