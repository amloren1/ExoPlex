
import numpy as np
import burnman
from scipy import interpolate

def get_rho(Planet,Mantle_filename,Core_wt_per,structural_params):

    radius_layers, rho_layers, rho0_layers, volume_layers, \
     mass_layers, cumulative_mass, Temperature_layers, \
     gravity, Pressure_layers, rhon, alpha, cp, phase=Planet

    num_mantle_layers, num_core_layers, number_h2o_layers = structural_params[4:]

    num_layers = num_mantle_layers+num_core_layers+number_h2o_layers

    for i in range(num_layers+1):
        if i < num_core_layers:
            rho_layers[i] = get_core_rho(Pressure_layers[i],Temperature_layers[i],Core_wt_per)
        elif i>=num_core_layers+num_mantle_layers:
            rho_layers[i] = get_mantle_rho(Pressure_layers[i],Temperature_layers[i],Mantle_filename)
        else:
            rho_layers[i] = get_water_rho(Pressure_layers[i],Temperature_layers[i])

    return rho_layers

def get_core_rho(Pressure,Temperature,Core_wt_per):
    frac_Si = Core_wt_per.get('Si')
    frac_O = Core_wt_per.get('O')
    frac_S = Core_wt_per.get('S')
    frac_Fe = Core_wt_per.get('Fe')


    return density

def get_mantle_rho(Pressure,Temperature,Mantle_filename):

    file = open(Mantle_filename+'_1.tab','r')
    temp_file = file.readlines()
    num_rows = len(temp_file[13:])
    num_columns = len(temp_file[12].split())


    data = temp_file[13:]
    grid = np.zeros((num_rows,num_columns))

    for i in range(num_rows):
        #for j in range(num_columns):
        columns = data[i].strip('\n').split()
        grid[i] = [float(j) for j in columns]


    pressure_grid = [row[1] for row in grid]
    temperature_grid = [row[0] for row in grid]
    density_grid = [row[2] for row in grid]

    density = interpolate.griddata((pressure_grid,temperature_grid),density_grid,(Pressure,Temperature),method='linear')
    return density


def get_water_rho(Pressure,Temperature):
    return 0