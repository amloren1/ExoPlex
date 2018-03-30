# coding: utf-8
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.


#####
import os
import sys
import numpy as np

import matplotlib.pyplot as plt
import pdb
#location from file which called this module
prem_dat = 'PREM.csv'

REarth = 6371 #kilometers

verbose = True

#Function just takes data from PREM
'''
To cite IRIS DMC Data Products effort:

Trabant, C., A. R. Hutko, M. Bahavar, R. Karstens, T. Ahern, and R. Aster (2012), 
 Data Products at the IRIS DMC: Stepping Stones for Research and Other Applications,
 Seismological Research Letters, 83(5), 846–854, doi:10.1785/0220120032.

To cite the source of this Earth model:

Dziewonski, A.M., and D.L. Anderson. 1981. “Preliminary reference Earth model.” Phys. Earth Plan. Int. 25:297-356.

To reference the use of this Earth model hosted by EMC:

IRIS DMC (2010), Data Services Products: PREM Preliminary Reference Earth Model, http://doi:10.17611/DP/9991844.
'''

def prem():
    #Defining each column
    #1 R,	Earth	radius	(km)
    #	2 Z,	Depth	(km)
    #	3 rho,	Density	(gm/cc)
    #	4 VPV,	P-wave	velocity	(km/sec)	{in	the	vertical	direction	where	anisotropic)
    #	5 VPH,	P-wave	velocity	(km/sec)	{in	the	horizontal	direction	where	anisotropic)
    #	6 VSV,	S-wave	velocity	(km/sec)	{in	the	vertical	direction	where	anisotropic)
    #	7 VSH,	S-wave	velocity	(km/sec)	{in	the	horizontal	direction	where	anisotropic)
    #	8 epsilon,	an	anisotropy	parameter	(ignore)
    #	9 Qmu	Quality	factor	relating	to	rigidity,	 (ignore)
    #	10 Qkappa.	Quality	factor	relating	to	incompressibility,	(ignore)
    
    cwd = os.path.dirname(__file__)

    dat = np.genfromtxt(cwd+'/PREM.csv', delimiter = ',', usecols = (0,1,2,3,4,5,6))
    
    #radius where r[0] is the center
    rad         = dat[:,0]
    depth       = dat[:,1]
    rho_depth   = dat[:,2]
    rho_rad     = rho_depth[::-1]

    
    prem_data = {'radius': rad, 'depth': depth, 'rho_depth': rho_depth, 'rho_radius': rho_rad , \
                   'VPV': dat[:,3], 'VPH': dat[:,4], 'VSV': dat[:,5], 'VSH': dat[:,6]}
    return(prem_data)
 

