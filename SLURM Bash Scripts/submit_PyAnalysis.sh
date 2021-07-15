#!/bin/bash -l
#SBATCH --job-name=PyAnalysis
# speficity number of nodes
#SBATCH -n 1

# specify the walltime e.g 20 mins
#SBATCH -t 10:00:00

# run from current directory
cd $SLURM_SUBMIT_DIR

module load anaconda

# command to use
# Write stdout+stderr to file
#SBATCH -o output.txt

ulimit -s unlimited

rm PythonLogFile

for a in $(seq 0.1 0.1 0.6) ; do
	echo "Parsing $a anisotropy data" >> PythonLogFile
	python Analysis.py -n 4 -o "HexOutput" -i "[20,40]" -a $a
done