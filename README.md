Square Lattice

The SimpleSquare.c can be compiled on Sonic using "gcc -o SimpleSquare SimpleSquare.c -lm".
The math.h library needs to be linked for successful compile.

Current Command Line Tags for SimpleSquare.c

-N : The dimension of the square lattice NxN. Default is 20.
-U : The upper bound on simulation steps.
-L : The lower bound on simulation steps.
-a : Set the Axis Anisotropy value. Default is 0.0
-e : Set the Exchange Anisotropy value. Default is 0.0
-i : Set the index for the output file. Needed creating multiple files. Default is 0.
-T : Set the estimate for the critical temperature. This allows for more sampled steps around the high variance region. Default is 1.2.
-s : Set Nearest Neighbours parameter. Default is 4.

INPUTVECS
Input file for SimpleSquare.c. This should list the nearest neighbours in terms of primitive vector coefficients around the centre point.
Example for square would be:
1,0
-1,0
0,1
0,-1


submit.sh
The submit.sh is an example file for submitting the SimpleSquare job to the Sonic cluster.

The SquareLatticeJob.sh file is for batching multiple ensembles for later averaging.


Analysis.py
This python script loads a batch of ensembles of the same dimension and creates the graph of average magnetism,
magnetic susceptibility and the U2 Binder Cumulant. This will be later converted into a Jupyter Notebook.


HeisenbergSquareLattice(...).py
These are depreciated python files for the SimpleSquare program. Will run but optimised and very slow.