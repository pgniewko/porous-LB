DESCRIPTION
==================================================
`porous-LB` is a code that performs Lattice-Boltzmann simulations for porous materials.

The code calculates material permeability (using [Darcy's law](https://en.wikipedia.org/wiki/Darcy%27s_law)) and [tortuosity](http://aip.scitation.org/doi/abs/10.1063/1.4711147) 
using BGK or MRT relaxation models.   
The calculations can be performed with periodic and non-periodic boundary conditions.   

GETTING THE CODE
==================================================
* To get the code:
```
git clone git@github.com:pgniewko/porous-LB.git
```

* To obtain the most recent version of the code:
```
git pull origin master
```

COMPILING AND INSTALLATION - UNIX
==================================================
* Executables
To build executables run:
```
make
```

Compilation options can be set by changing proper variables in Makefile:

* `palabosRoot` points to the directory with Palabos source code   
* `compileFlags` sets the compilation flags. Uncomment `-DMRT` option for MRT dynamics


EXTERNAL LIBRARIES
================
The code requires an external library: [Palabos](http://www.palabos.org/).

USAGE
=====
Run the code with no parameters for the options list and usage:
```
./permeability 

```
The arguments are:
    
1.  Input file name
2.  Output directory name
3.  Output file name
4.  Number of lattice sites in X direction
5.  Number of lattice sites in Y direction
6.  Number of lattice sites in Z direction
7.  Pressure difference delta P
8.  Direction of the pressure gradient
9.  Periodicity flag. 
    * 1 for PBC in a direction perpendicular to argument no.8
10. Refinemenet level - needed for the max iterations estimation

INPUT FILE
==========
The input file (the 1st argument) contains a single line.   
The file encodes solid and liquid phase grids in a binary manner: 1 for solid, and 0 for liquid phases.    
The grid points are looped over X,Y, and Z directions.   

For example, for 3x3x3 cube filled with liquid(0), and a solid(1) site in the middle of this cube, we can create an inpute file with the following, simple script in Python:    
```Python
import numpy as np

lattice = np.zeros([3,3,3])
lattice[1][1][1] = 1
s = ""

for i in range(3):
  for j in range(3):
    for k in range(3):
      s += str( lattice[i][j][k]) + " "

print(s)
```

COPYRIGHT NOTICE
================
Copyright (C) 2017-2018, Pawel Gniewek   
Email : gniewko.pablo@gmail.com   
All rights reserved.   
License: BSD 3   

ACKNOWLEDGMENTS
===============

The code is developed based on the Palabos' on-line [example](http://www.palabos.org/documentation/tutorial/permeability.html).


