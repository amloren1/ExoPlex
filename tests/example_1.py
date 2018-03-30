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

import pdb

print sys.path

import exoplex.run as run




#####
'''
testing grids functions
'''

####

FeMg = np.arange(0, .2,0.1)
SiMg = np.arange(0.1,.3,0.1)
CMF  = np.arange(0.3,0.5,0.1)

run.grid_cmf(mass = 1.0,femg = FeMg, simg = SiMg, cmf = CMF, filename = 'testing1212.dat')

sys.exit()

#Use inputs python file to call exoplex
# this creates an array of planet model data that you may use
# to plot or output files
#**NOTE: ENTER FILENAME WITHOUT .py**
# inputs_1 asks exoplex to model one planet. 

Planets = run.exoplex('input_1')
pdb.set_trace()

#all model data is now in Planets. Lets make some plots and print the 
# results to a file

run.write(planet = Planets, filenames = ['planet_1.dat'])


#our first output will be a plot against the PREM (Dziewonski & Anderson 1981)

#run.pltprem(planet = Planets,label = ['planet 1', 'planet 2'])

#run.pltrho(planet = Planets,label = ['planet 1', 'planet 2'])


