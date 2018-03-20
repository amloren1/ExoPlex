# coding: utf-8
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.


#**********************************************************************#
'''
'''
#**********************************************************************#



import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))
    
import pdb
import run 



#####
'''
create a data grid for a planet of some mass
vary mantle Fe/Mg, Si/Mg, water mass fraction, and the CMF (core mass fraction)
This is good to make if there is not a premade grid for this 
specific range of mantle composition
or core mass fraction. 

Output will appear in Solutions/Grids folder. 

kwargs:
mass     = x  #mass of planet in Earth masses
femg     = [star, stop, step] range of elemental ratio of Fe to Mg in the mantle
simg     = [star, stop, step] range of elemental ratio of Si to Mg in the mantle
cmf      = [star, stop, step] core mass fraction range
h2o      = [star, stop, step] (optional) range of water mass fraction
filename = x                  (optional) name of file to be stored in Solutions/Grids 
'''
####

FeMg = [0.1, 0.3,0.1]
SiMg = [0.2, 0.4,0.1]
CMF  = [0.3, 0.5,0.1]
wt_h2o = [0.0, 0.0, 0.1]
run.cmf_grid(mass = 1.0,femg = FeMg, simg = SiMg, cmf = CMF, filename = 'testing1212.dat')

sys.exit()











