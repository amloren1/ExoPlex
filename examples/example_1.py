# coding: utf-8
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.


#**********************************************************************#
'''
this example script will show how to call exoplex and create models
based off of the inputs from some input file (open inputs_1.py for corresponding
file). 


***Before you run this script***
-use pip to install exoplex from PYPI
-create a directory structure for exoplex to work with
	-this script should sit alongside a directory called Solutions/
	and Solutions/Grids
	-you need to have params.py and some input file (use input_1.py as a template) next to
	this script
	-example directory tree if you have this script in a directory called ExoPlex/

	ExoPlex/
		-example_1.py
		-params.py
		-input_1.py
		Solutions/
			#mantle thermodynamic data for each composition
			Grids/
				#grid results will be placed here

'''
#**********************************************************************#



import sys

import exoplex.run as run




#####
'''
function: exoplex()

Use inputs from python file to call exoplex

returns:
arrays of model data. Each element in the array will correspond to a dictionary of data for
that model. 

arg:
-name of an input python file
	- must be in the same format as input_1.py
	- just type the name of the file without the .py
	- input file must be stored beside the script you call
	exoplex from

'''
#####


#**NOTE: ENTER FILENAME WITHOUT .py**

Planets = run.exoplex('input_1')




#####
'''
function: write(**kwargs)

Writes model results to a file. 

returns:
file with results

kwargs:
planet    = dictionary array returned from run.exoplex
filenames = (optional) list of strings that represent desired filenames for each planet model
			defaults to planet0*.dat
'''
#####


run.write(planet = Planets, filenames = ['planet_1.dat'])


#####
'''
function: pltprem(**kwargs)

Plot exoplex models against the PREM (Dziewonski & Anderson 1981) in 
a depth vs. density diagram

returns:
file with results

kwargs:
planet    = dictionary array returned from run.exoplex
label     = (optional) list of strings for labels on the plot
'''
#####



run.pltprem(planet = Planets,label = ['planet 1'])


#####
'''
function: pltrho(**kwargs)

Plot exoplex models of depts vs density

returns:
file with results

kwargs:
planet    = dictionary array returned from run.exoplex
label     = (optional) list of strings for labels on the plot
'''
#####


run.pltrho(planet = Planets,label = ['planet 1'])


