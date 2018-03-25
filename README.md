
### run.mvr_grid_cmf(**kwargs)
Solves for the mass and radius relationship for a range of compositions and masses. Solutions are output as data files and a mass vs. radius plot. Each combination of composition is given a data file. the cmf designation signifies that this function takes in mantle composition and core mass fraction as opposed to mvr_grid, which calculates core mass fraction from bulk elemental abundances.  

##### kwargs
 * mass    : python list of mass range in the format [start, stop, step]  
 * femg    : python list of bulk (global) elemental Fe to Mg ratio range in the format [start, stop, step] 
 * simg    : python list of bulk (global) elemental Si to Mg ratio range in the format [start, stop, step] 
 * cmf     : python list of core mass fractions in the format [start, stop, step]
 * plot    : True or False. Plots mass vs radius for composition range if True
 * h2o     : (optional) range of water mass fraction in the format [start,stop,step]
 * filenames: (optional) list of filenames for each composition. Otherwise default filenames will be created
 
 
##### returns
 
 * multiple files with solutions for a range of masses. Each composition creates a file. 
 
##### example

```python
M      = [0.3, 2.5, 0.1]
FeMg   = [0.9, 1.5, 0.1]
SiMg   = [1.1, 1.5, 0.1]
CMF    = [0, 0.6, 0.05]
wt_h2o = [0, 0.2, 0.05]

run.mvr_grid_cmf(femg = FeMg, simg = SiMg, cmf = CMF, h2o = wt_h2o)
```
###### returns
![alt text](https://i.imgur.com/OgyabcA.png)
