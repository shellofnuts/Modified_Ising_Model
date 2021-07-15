# Magnetic Modelling using Heisenberg Model

This code implements a monte-carlo approach to simulating magnetic moments in 2D lattices. The N=3 Heisenberg model is implemented with anisotropy factors.
This allows us to overcome the Mermin-Wagner theorem and stimulate and tune long-range magnetic order in 2D lattices.
The lattice symmetry can be arbitrarily defined, where the strength and interaction type (FM/AFM) can be individually defined for each direction.

## Compiling

The LatticeSim_OpenMPI.c can be compiled using:
```
gcc foo.c -o foo -lm
```

The math.h library needs to be linked for successful compile.

## Current Command Line Tags

These

- -N : REQUIRED: The dimension of the square lattice NxN. Default is 20.
- -U : REQUIRED: The upper bound on simulation steps.
- -L : DEPRECIATED: The lower bound on simulation steps.
- -a : Set the Axis Anisotropy value. Default is 0.0
- -e : Set the Exchange Anisotropy value. Default is 0.0
- -i : Set the index for the output file. Needed for creating multiple files. Default is 0.
- -T : Set the estimate for the critical temperature. This allows for more sampled steps around the high variance region. Default is 1.2.
- -s : DEPRECIATED: Set Nearest Neighbours parameter. Default is 4. This should be set in the INPUTVECS file
- -I : Set the interval between each temperature point.
- -r : Set the temperature range of the simulation.



## INPUTVECS
Input file for simulation. This should list the nearest neighbours in terms of primitive vector coefficients around the centre point.

The format of the file is given as:
1. Simulation name
2. Number of basis sets, Number of nearest neighbours
3. Basis set number
4. X, Y, Interaction Type
5. Repeat 4 for all nearest neighbours
6. Repeat 3-5 for all basis sets

An exampled for a hexagonal lattice is given below:

```
Hexagonal Lattice
2,3
1
0,0,1
0,1,1
-1,1,1
2
0,0,1
0,1,1
-1,1,1
```

## submit.sh
The submit.sh is an example file for submitting the SimpleSquare job to the Sonic cluster.

The SquareLatticeJob.sh file is for batching multiple ensembles for later averaging.


Analysis.py
This script extracts and averages 3 data points of interest: the average magnetisation across multiple ensemble averages, the average magnetic susceptibility, and the average U2 Binder Cumulant.

File structure should have a directory with folders titiled {}x{}_Data where {} is the lattice size and contains the ensemble averages at that lattice size.

Commmand Line arguments are:
-f: Root of file structure discussed above.
-n: Number of .csv files to average over.
-o: The output file name.
-i: The lattice sizes to loop over. Must be passed as a list enclosed in "". e.g: "[20, 40, 60]".


Lattice_Analysis.ipynb
This JupyterNotebook reads in the output CSV file and plots them on inified graphs.
