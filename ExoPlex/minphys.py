
import numpy as np
import burnman
from scipy import interpolate

def get_rho(Planet,mass_planet,grids,Core_wt_per,structural_params):

    radius_layers, rho_layers, rho0_layers, volume_layers, \
     mass_layers, cumulative_mass, Temperature_layers, \
     gravity, Pressure_layers, rhon, alpha, cp, phase=Planet

    num_mantle_layers, num_core_layers, number_h2o_layers = structural_params[4:]

    num_layers = num_mantle_layers+num_core_layers+number_h2o_layers
    temperature_grid, pressure_grid, density_grid, speed_grid, alpha_grid, cp_grid, phases_grid = grids


    for i in range(num_layers+1):
        if i < num_core_layers:
            rho_layers[i] = get_core_rho(Pressure_layers[i],Temperature_layers[i],Core_wt_per)
        elif i>=num_core_layers+num_mantle_layers:
            rho_layers[i] = get_mantle_rho(Pressure_layers[i],Temperature_layers[i],pressure_grid,\
                                           temperature_grid,density_grid)
        else:
            rho_layers[i] = get_water_rho(Pressure_layers[i],Temperature_layers[i])

    return rho_layers

def get_core_rho(Pressure,Temperature,Core_wt_per):
    frac_Si = Core_wt_per.get('Si')
    frac_O = Core_wt_per.get('O')
    frac_S = Core_wt_per.get('S')
    frac_Fe = Core_wt_per.get('Fe')


    return density

def get_mantle_rho(Pressure,Temperature,pressure_grid,temperature_grid,density_grid):
    density = interpolate.griddata((pressure_grid,temperature_grid),density_grid,(Pressure,Temperature),method='linear')
    return density


def get_water_rho(Pressure,Temperature):
    return 0