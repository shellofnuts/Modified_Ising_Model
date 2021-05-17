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
This script extracts and averages 3 data points of interest: the average magnetisation across multiple ensemble averages, the average magnetic susceptibility, and the average U2 Binder Cumulant.

File structure should have a directory with folders titiled {}x{}_Data where {} is the lattice size and contains the ensemble averages at that lattice size.

Commmand Line arguments are:
-f: Root of file structure discussed above.
-n: Number of .csv files to average over.
-o: The output file name.
-i: The lattice sizes to loop over. Must be passed as a list enclosed in "". e.g: "[20, 40, 60]".


Lattice_Analysis.ipynb
This JupyterNotebook reads in the output CSV file and plots them on inified graphs.
