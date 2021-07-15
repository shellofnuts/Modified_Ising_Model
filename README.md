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

## Python Analysis

### Analysis.py

The Analysis.py requires Pandas, Numpy & Scipy.stats.moment() to parse the .csv output files.
It reads in the raw .csv files, finds the average magnetic moment and standard deviation at each temperature point.
The magnetic susceptibility and U2 Binder cumulant for the lattice size is created.
The output file is of the following format.

| Temperature | {}x{} Average Magnetic Moment | {}x{} Average Susc | {}x{} Average U2 | {}x{} Magnetic Moment std | Repeated for larger cell sizes |
| ----------- | ----------------------------- | ------------------ | ---------------- | ------------------------- | ------------------------------ |

{} denotes the lattice size. The lattice sizes to work over are defined using the command line arguments:

- -f : Path to directory containing the input files.
- -n : Number of .csv files to average over for each setting value.
- -o : The output file name.
- -i : The lattice sizes to loop over. Must be passed as a list enclosed in "". e.g: "[20, 40, 60]".
- -a : Number denoting the file number. Can be used to express the used anisotropy value.

Input file format: ```python '{}_{}_spinDist_{}_{}.csv'.format(lattice_size, lattice_size, a, n) ```
Output file format: "-o"_"-a".csv


### Lattice_Analysis.ipynb

This Jupyter Notebook reads in the output .csv file, unifies the different lattices sizes onto three different graphs.

![Magnetic Transition Curve](/Results/MagneticCurve.png)

![Magnetic Susceptibility](/Results/MagneticSusceptibility.png)

![U2 Binder Cumulant](/Results/MagneticBinderCumulant.png)

It automatically formats the graphs to be outputted.


## Submitting to SLURM

The submit_Combined.sh is an example file for submitting the compiled job to the UCD Sonic cluster.
It is also an example of how to write a bash script to create a large number of ensembles and then create the output/averaged script.

The separate submit_PyAnalysis.sh and submit_MPI.sh are separate scripts to run the Python Analysis.py script and the LatticeSim_OpenMPI script respectively.

The headers of these scripts are particular to the UCD Sonic Cluster and its implementation of SLURM.
