

#**********************************************************************#
'''
This script enables the user to input a simplified composition 
Fe/Mg and Si/Mg, and a mass to find the corresponding radius, CMF, CRF
'''
#*********************************************************************



import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))


import ExoPlex as exo
import PREM.prem as p
from params import *
import pdb



################
'''
Steps:
1. import data: each file is for a mass which, in itself if the title name

2. store masses into a 1-D array. This will be the interpolation array

3. create routine to find the radii which correspond to input composition
'''
################
def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def bissect(grid, value):
    
    mid = int(len(grid)/2.)
    hi  = int(len(grid)-1)
    lo  = 0
    if value/max(grid) > 1.0+1e-6 :
        
        print 'Exiting. The following value is not in the grids yet:'
        print '{}'.format(value)
        print 'Max = {}'.format(max(grid))
        print 'ratio = {}'.format(value/max(grid))
        sys.exit()
    for i in range(len(grid)):
        if value > grid[mid]:
            lo = mid
            mid = int((hi+lo)/2.)
            #print mid
        elif value < grid[mid]:
            hi  = mid
            mid = int((hi+lo)/2.)
            #print mid
    return mid, lo, hi

def interpolant(lo, hi, val):
    q = (val-lo+0.0)/(hi-lo)
    return q

def interpd_data(q, y1, y2):
    
    y = y1+q*(y2-y1)
    
    return y




#Find row which corresponds to composition
# Fe/Mg is from 0.7-1.6, Si/Mg is 0.7-1.6

FeMg = np.arange(0.7, 1.7, 0.1)
SiMg = np.arange(0.7, 1.7, 0.1)

def loc_data(FeMg_grid, SiMg_grid, femg, simg):
    
    femg_i = bissect(FeMg_grid, femg)[0]
    simg_i = bissect(SiMg_grid, simg)[0]
    #loc is the row in any mass file which is the specified composition
    loc    = femg_i*len(FeMg_grid)+simg_i

    return loc
###
#File names for the data files
###
def file_names(mass):
    mass_string = repr(round(mass,2)).replace('.', ',')
    name = mass_string+'.dat'
    
    return name

########################################################################
'''
Import all of the mass files
'''
########################################################################

#mass grid for interpolation
Mass = np.arange(0.5, 1.8,0.1)

Rads = np.zeros(len(Mass))


####
#Pull data from arrays for a specific composition
####
def make_data_arrays(femg, simg):

    N_mas = len(Mass)
    t_femg = np.zeros(N_mas)
    t_simg = np.zeros(N_mas)
    rho    = np.zeros(N_mas)
    cmf    = np.zeros(N_mas)
    crf    = np.zeros(N_mas)
    #find row for this composition
    row = loc_data(FeMg, SiMg, femg, simg)
    #pull data row from each file
    for i in range(len(Mass)):
        f_name = file_names(Mass[i])
        #'Mass', 'Radius','Fe/Mg','Si/Mg', 'Bulk density', 'core_wt%', 'core_rad%'
        dat = np.genfromtxt(f_name)
    
        dat_row = dat[row,:]
    
        Rads[i]    = dat_row[1] 
        t_femg[i]  = dat_row[2]
        t_simg[i]  = dat_row[3]
        rho[i]     = dat_row[4]
        cmf[i]     = dat_row[5]
        crf[i]     = dat_row[6]
        
    keys = ['Mass', 'Radius','FeMg','SiMg', 'rho', 'CMF', 'CRF']
    
    return(dict(zip(keys,[Mass, Rads, t_femg, t_simg, rho, cmf, crf])))
   
####
# plot limits: pure fe, pure mantle and pure water/ice
#### 

def plot_limits(f_name):
    
    dat = np.genfromtxt(f_name)
    mas = dat[:,0]
    rads = dat[:,1]

    ax.plot(mas, rads, lw = 5, color = 'dimgray', alpha = 0.8)   
    ax.annotate('100% Fe', xy=(np.median(mas), np.median(rads)*0.978), \
     xytext = (np.median(mas), np.median(rads)*0.978), fontsize = 16, \
     rotation=10, weight = 'bold')

####
# take composition, find data arrays, interpolate to find radius
# plot on MvR
#### 
def plot_comps(femg, simg, mass, err):
    import random
    r = lambda: random.randint(0,255)

    #get data arrays
    def find_radius(femg, simg, mass):
        dat = make_data_arrays(femg, simg)
        Rads = dat.get('Radius')
        rho = dat.get('rho')
        cmf = dat.get('cmf')
        crf = dat.get('crf')
        
        #find interpolation value
        mid, lo, hi = bissect(Mass, mass)
        
        q = interpolant(Mass[lo],Mass[hi],mass)
     
        radius = interpd_data(q, Rads[lo], Rads[hi])
        
        return radius, Rads
    
    
    plot_limits('pure_Fe.dat')
    
    lab_size = 23
    tic_size = 18

    ax.set_xlabel(r"Mass (M$_\oplus$)", fontsize = lab_size )
    ax.set_ylabel(r"Radius (R$_\oplus$)", fontsize = lab_size)
    ax.tick_params(direction='in', length=6, labelsize = tic_size)
    ax.grid(color='grey', linestyle='-', alpha = 0.4, linewidth=.7)
    ax.set_xlim(min(Mass)*.9, max(Mass)*1.06)
    ax.set_ylim(0.5, 1.28)
    for i in range(len(femg)):
        
        radius, rads = find_radius(femg[i], simg[i], mass[i])
        no_lab = "_nolegend_"
        l1 = 'Fe/Mg = {}; Si/Mg = {}'.format(femg[i],simg[i])
        if i < len(femg)-1:
            l2 = 'Fe/Mg = {}; Si/Mg = {}'.format(femg[i+1],simg[i+1])
        else:
            l2 = 'dummd'
        
        j = ('#%02X%02X%02X' % (r(),r(),r()))
        ax.scatter(mass[i],radius, marker = 'H', linewidth = 7, c = j, edgecolor=j,  \
                  label = 'R(M={})={}'.format(mass[i], radius))
        ax.errorbar(mass[i], radius, yerr=0,xerr = err[i], \
            fmt='o', elinewidth=3, capthick=3, color = j)


        if l1 != l2:
            ax.plot(Mass, rads, label = l1, lw = 5, alpha = 0.8)
        plt.draw()
    
    
    
    plt.legend(loc = 'lower right', fontsize = tic_size, scatterpoints=1)
    plt.show()

    return radius


fig, ax = plt.subplots(figsize = (20,10))
plt.rc('font', family='serif')
plot_comps([0.7,1.4, 1.2,1.2], [1.2,1.6, 1.5, 1.5], [0.67, 1.2, 0.5, 1.0], \
               [0.07, 0.2, 0.09, 0.5])


pdb.set_trace()
sys.exit()


test = 1.647
mid, lo, high = bissect(Mass, test)
Q = interpolant(Mass[lo], Mass[high], test)
print '\n{} is at location = {} '.format(test, mid)
print 'mid = {} \nlow = {} \nhi = {}\n'.format(mid, lo, high)
print 'interpolant between M_lo = {} and M_hi = {} \nis q = {}'.format(Mass[lo], \
                                Mass[high],Q) 

print Mass  










