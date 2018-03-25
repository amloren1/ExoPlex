# About ExoPlex
A self-consistent mineral physics solver for exoplanets. We model rocky exoplanets with a core, lower and upper mantle, and water/ice layer. We use the Gibbs free energy minimization package, Perple_x (Connolly 09) along with a solution method from Stixrude 2011 to model mantles with natural phase transformations at depth (Pressure and Temperature dependence). 

#### Current authors:
- Cayman Unterborn
- Alejandro Lorenzo

#### Science contributors:
- Steven Desch
- S.-H. Dan Shim: mineral physics
- Byeongknwan Ko: mineral physics

## Getting Started

ExoPlex is written in Python 2.7 and is currently stable to its full capability in the Mac OS and Linux platforms. The mineral physics solver for mantle phase transformations relies on a Fortran program, perple_x which is freely available from their website. We have included some precompiled phase solutions from perple_X so exoplex can be run without installing perple_x however, the user is limited in choice of planet composition to what is already included. Below we detail the full dependencies for exoplex.

### Requirements 

- Mac OS or Linux/Unix

- Perplex_x version 6.8.1 (updated 3.16.18)

- Python 2.7+

- Python libraries (use latest versions):
  - numpy: handle arrays 
  - os: search directories for files
  - math: math package 
  - matplotlib: plotting results from ExoPlex
  - pexpect: used to call Perple_X from python
  - multiprocessing: not absolutely necessary but if you can get this working, it cuts code runtime in half 
  - scipy: interpolation and integration
  - burnman: solves equations of state BME II,III, IV for the core and ice layers

 
#### Installing on Ubuntu and its derivatives

1. Clone the master branch of this repository. 
2. Extract tar file within the ```Solutions/``` directory:

```
tar -xvzf solutions.tar.gz
```

3. In a directory above ExoPlex, or otherwise not within the ExoPlex directory, install the latest version of Perple_X
    - [Perple_X linux version](http://www.perplex.ethz.ch/perplex/ibm_and_mac_archives/LINUX/Links_for_the_latest_LINUX_version_of_Perple_X.html) 
    


4. Move the following executable programs into the ```ExoPlex/PerPlex/``` directory:
   - build, vertex, werami
   - the other files from the Perple_x install will not be used by ExoPlex

5. Install the libraries mentioned above. All of these are available in the PyPI and can be installed with pip, 
```
pip install package_name
```


#### Installing on MAC OS

1. Clone the master branch of this repository. 
2. Extract tar file within the ```Solutions/``` directory:

```
tar -xvzf solutions.tar.gz
```


3. In a directory above ExoPlex, or otherwise not within the ExoPlex directory, install Perple_X version 6.8.1

    - [Perple_X MAC OS version](http://www.perplex.ethz.ch/perplex/ibm_and_mac_archives/OSX/) 
    

4. Move the following executable programs into the ```ExoPlex/PerPlex/``` directory:
   - build, vertex, werami
   - the other files from the Perple_x install will not be used by ExoPlex

5. Install the libraries mentioned above. All of these are available in the PyPI and can be installed with pip, 
```
pip install package_name
```
## Using ExoPlex

ExoPlex is meant to be run as a library complimentary to python scripts. The idea is to store your scripts in a directory on the same level as the ```ExoPlex/``` directory. Most of the functions the user will be interacting with are located in the ```run/``` directory. To use these from the  ```start_here/``` directory, for example, scripts should begin with:

```python
import os
import sys
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))
import run
```

Where the third line alters the path to enable the script to run functions from the other directories. 





## Functions

### run.exoplex(inputs_file)
Calls ExoPlex functions and executes desired models based on input python script. 

##### args
 * inputs_file: String, name of a python script with input parameters. Example format is given in ```start_here/input_example.py```
 
##### returns
 * numpy array of dictionaries which contain model information with keys: 
 ['volume', 'phases', 'phase_names', 'temperature', 'density', 'dmass', 'K', 'mantle_ratios', 'gravity', 'Vp', 'pressure', 'Vs', 'radius', 'Vphi', 'bulk_ratios', 'alpha', 'cp', 'mass']

 
##### example

```python
planets = run.exoplex('inputs_1')
```





### run.pltprem(**kwargs)
Radial density profile of each modeled planet plotted with the PREM (A. M. Dziewonski & D. L. Anderson 1981).

##### kwargs
 * planet: results from `run.exoplex()`  
 * label: string array (or list), label for each model to be plotted. This will be displayed in a legend on the plot. 
 
##### returns
 
 * None
 
##### example

```python
planets = run.exoplex('inputs_1')
run.pltrho(planet = planets, label = ['Mg/Si=1', 'Mg/Si=2'])
```
 
###### returns
![alt text](https://i.imgur.com/YvZ1FT8.png)







### run.pltrho(**kwargs)
Radial density profile of each modeled planet.

##### kwargs
 * planet: results from `run.exoplex()`  
 * label: string array (or list), label for each model to be plotted. This will be displayed in a legend on theplot. 
 
##### returns
 
 * plot
 
##### example

```python
planets = run.exoplex('inputs_1')
run.pltrho(planet = planets, label = ['Mg/Si=1', 'Mg/Si=2'])
```
###### returns
![alt text](https://i.imgur.com/rPfPOGG.png)


### run.single_grid(**kwargs)
Creates a data file for a planet of some mass for a range of bulk compositions. Core mass fraction is calculated from the inputs. 

##### kwargs
 * mass: mass in units of Earth mass  
 * femg: python list of bulk (global) elemental Fe to Mg ratio range in the format [start, stop, step] 
 * simg: python list of bulk (global) elemental Si to Mg ratio range in the format [start, stop, step] 
 * h2o: (optional) range of water mass fraction in the format [start,stop,step]
 * filename: (optional) name of file, otherwise default filename will be created
 
 
##### returns
 
 * file of solutions for radius as a function of mass for the desired composition range
 
##### example

```python
FeMg   = [0.9, 1.5, 0.1]
SiMg   = [1.1, 1.5, 0.1]
wt_h2o = [0, 0.2, 0.05]

run.single_grid(femg = FeMg, simg = SiMg, h2o = wt_h2o, filename = 'composition_1.dat')
```




### run.single_grid_cmf(**kwargs)
Creates a data file for a planet of some mass for a range of compositions. Similar to single_grid however the user specifies the mantle composition directly along with the core mass fraction. 

##### kwargs
 * mass    : mass in units of Earth mass  
 * femg    : python list of bulk (global) elemental Fe to Mg ratio range in the format [start, stop, step] 
 * simg    : python list of bulk (global) elemental Si to Mg ratio range in the format [start, stop, step] 
 * cmf     : python list of core mass fractions in the format [start, stop, step]
 * h2o     : (optional) range of water mass fraction in the format [start,stop,step]
 * filename: (optional) name of file, otherwise default filename will be created
 
 
##### returns
 
 * file of solutions for radius as a function of mass for the desired composition range
 
##### example

```python
FeMg   = [0.9, 1.5, 0.1]
SiMg   = [1.1, 1.5, 0.1]
CMF    = [0, 0.6, 0.05]
wt_h2o = [0, 0.2, 0.05]

run.single_grid(femg = FeMg, simg = SiMg, cmf = CMF, h2o = wt_h2o, filename = 'composition_1.dat')
```




### run.mvr_grid(**kwargs)
Solves for the mass and radius relationship for a range of compositions and masses. Solutions are output as data files and a mass vs. radius plot. Each combination of composition is given a data file. The core mass fraction is calculated from the bulk elemental abundances of the planet. 

##### kwargs
 * mass    : python list of mass range in the format [start, stop, step]  
 * femg    : python list of bulk (global) elemental Fe to Mg ratio range in the format [start, stop, step] 
 * simg    : python list of bulk (global) elemental Si to Mg ratio range in the format [start, stop, step] 
 * plot    : True or False. Plots mass vs radius for composition range if True
 * h2o     : (optional) range of water mass fraction in the format [start,stop,step]
 * filenames: (optional) list of filenames for each composition. Otherwise default filenames will be created
 
 
##### returns
 
 * multiple files with solutions for a range of masses. Each composition creates a file. 
 * a plot of mass vs. radius
 
##### example

```python
M      = [0.3, 2.5, 0.1]
FeMg   = [0.9, 1.5, 0.1]
SiMg   = [1.1, 1.5, 0.1]
wt_h2o = [0, 0.2, 0.05]

run.mvr_grid(femg = FeMg, simg = SiMg, h2o = wt_h2o, plot =False)
```





### run.mvr_grid_cmf(**kwargs)
Solves for the mass and radius relationship for a range of compositions and masses. Solutions are output as data files and a mass vs. radius plot. Each combination of composition is given a data file. the cmf designation signifies that this function takes in mantle composition and core mass fraction as opposed to mvr_grid, which calculates core mass fraction from bulk elemental abundances.  

##### kwargs
 * mass    : python list of mass range in the format [start, stop, step]  
 * femg    : python list of mantle elemental Fe to Mg ratio range in the format [start, stop, step] 
 * simg    : python list of mantle elemental Si to Mg ratio range in the format [start, stop, step] 
 * cmf     : python list of core mass fractions in the format [start, stop, step]
 * plot    : True or False. Plots mass vs radius for composition range if True
 * h2o     : (optional) range of water mass fraction in the format [start,stop,step]
 * filenames: (optional) list of filenames for each composition. Otherwise default filenames will be created
 
 
##### returns
 
 * multiple files with solutions for a range of masses. Each composition creates a file. 
 * a plot of mass vs. radius 
 
##### example

```python
M      = [0.3, 2.5, 0.1]
FeMg   = [0.2, 0.3, 0.1]
SiMg   = [0.3, 0.5, 0.1]
CMF    = [0.3, 0.5, 0.1]
wt_h2o = [0, 0.1, 0.1]

run.mvr_grid_cmf(femg = FeMg, simg = SiMg, cmf = CMF, h2o = wt_h2o, plot = True)
```
###### returns
![alt text](https://i.imgur.com/OgyabcA.png)






### run.write(**kwargs)
Radial density profile of each modeled planet.

##### kwargs
 * planet: results from `run.exoplex()`  
 * filenames: string array (or list). Filename for each planet modeled.  
 
##### returns
 
 * Prints to file a grid of mass, radius, density, pressure, and temperature, at each shell. Mass and radius are accumulative values from the planet center outwards. 
 
##### example

```python
planets = run.exoplex('inputs_1')
run.write(planet = planets, filenames = ['planet_1.dat', 'planet_2.dat'])
```






### run.writeall(**kwargs)
Radial density profile of each modeled planet.

##### kwargs
 * planet: results from `run.exoplex()`  
 * filenames: string array (or list). Filename for each planet modeled.  
 
##### returns
 
 * Prints to file a grid of mass, radius, density, pressure, temperature, heat capacity, thermal emissivity, and mineral profile of the compound at each shell. Mass and radius are accumulative values from the planet center outwards. 
 
##### example

```python
planets = run.exoplex('inputs_1')
run.writeall(planet = planets, filenames = ['planet_1.dat', 'planet_2.dat'])
```




## Examples

Learning how to use ExoPlex is best done by working through the examples provided in the ```start_here/``` directory. A new user should be able to gain enough insight from the current examples to use ExoPlex to its fullest however, we will be adding more exampes in the future.


## Contributing

ExoPlex is open source and we encourage users to contribute. In this case, fork the project and make a pull request. 



## License

Copyright (C) 2017 - by the ExoPlex team, released under the GNU
GPL v2 or later.

## Acknowledgments

* James Connolly (PerPle_X) for providing tips and source code
* This project was done as part of the NASA NExSS grant

## References 

Connolly JAD (2009) The geodynamic equation of state: what and how. Geochemistry, Geophysics, Geosystems 10:Q10014 DOI:10.1029/2009GC002540.

