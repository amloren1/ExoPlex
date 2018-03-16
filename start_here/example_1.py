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


#Use inputs python file to call exoplex
# this creates an array of planet model data that you may use
# to plot or output files
#**NOTE: ENTER FILENAME WITHOUT .py**
# inputs_1 asks exoplex to model one planet. 

Planets = run.exoplex('input_2')


#all model data is now in Planets. Lets make some plots and print the 
# results to a file

run.write(planet = Planets, filenames = ['planet_1.dat', 'planet_2.dat'])


#our first output will be a plot against the PREM (Dziewonski & Anderson 1981)
#run.pltprem(planet = Planets)

run.pltrho(planet = Planets,label = ['planet 1', 'planet 2'])


