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

ExoPlex is written in Python 2.7 and is currently stable to its full capability in the Mac OS and Linux platforms. The mineral physics solver for mantle phase transformations relies on a Fortran program, Perple_x which is freely available from their website.
We have implemented a version with precompiled phase diagram solutions and results from ExoPLex which can be searched on a grid. This version is located in another repository called ExoPlex Gridsearch. This version is platform independent.
Below we detail the dependencies for the full version of ExoPLex which does rely on Perple_x. Using precompiled phase solutions does not depend on Perple_X so the user may skip those steps if they wish.

### Requirements 

- Mac OS or Linux/Unix

- Perplex_x version 6.7

- Python 2.7+

- Python libraries (use latest versions):
  - numpy: handle arrays 
  - os: search directories for files
  - math: math package 
  - matplotlib: plotting results from ExoPlex
  - pexpect: used to call Perple_X from python
  - multiprocessing: not absolutely neccesary but if you can get this working, it cuts code runtime in half 
  - scipy: interpolation and integration
  - burnman: solves equations of state BME II,III, IV for the core and ice layers

 
### Installing on Ubuntu and its derivatives

1. Clone the master branch of this repository. 
2. Extract tar file within the ```Solutions/``` directory:

```
tar -xvzf solutions.tar.gz
```

3. In a directory above ExoPlex, or otherwise not within the ExoPlex directory, install the latest version of Perple_X
    - [Perple_X linux version](http://www.perplex.ethz.ch/perplex/ibm_and_mac_archives/LINUX/Links_for_the_latest_LINUX_version_of_Perple_X.html) 
    


4. Move the following executable programs into the ```ExoPlex/PerPlex/``` directory:
   - build, vertex, werami
   - e.g. from within the ExoPlex/PerPlex   ```path/to/perplex/install .```
   - the other files from the Perple_x install will not be used by ExoPlex

5. Install the libraries mentioned above. All of these are available in the PyPI and can be installed with pip, ```pip install package```


### Installing on MAC OS

1. Clone the master branch of this repository. 
2. Extract tar file within the ```Solutions/``` directory:

```
tar -xvzf solutions.tar.gz
```

3. In a directory above ExoPlex, or otherwise not within the ExoPlex directory, install the latest version of Perple_X
    - [Perple_X MAC OS version](http://www.perplex.ethz.ch/perplex/ibm_and_mac_archives/OSX/) 
    

4. Move the following executable programs into the ```ExoPlex/PerPlex/``` directory:
   - build, vertex, werami
   - e.g. from within the ExoPlex/PerPlex   ```path/to/perplex/install .```
   - the other files from the Perple_x install will not be used by ExoPlex

5. Install the libraries mentioned above. All of these are available in the PyPI and can be installed with pip, ```pip install package```


## Examples

Learning how to use ExoPlex is best done by working through the examples provided in the ```Examples/``` directory. A new user should be abe to gain enough insight from the current examples to use ExoPlex to its fullest however, we will be adding more exampes in the future.

### Earth_models.py

### autonio.py

## Contributing

ExoPlex is open source and we encourage users to contribute. In this case, fork the project and make a pull request. 



## License

Copyright (C) 2017 - by the ExoPlex team, released under the GNU
GPL v2 or later.

## Acknowledgments

* James Connolly (PerPle_X) for providing tips and source code
* This project was done as part of the NASA NExSS grant
