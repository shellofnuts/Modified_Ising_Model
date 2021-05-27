#!/bin/bash -l
#SBATCH --job-name=CSQ_Spin1
# speficity number of nodes 
#SBATCH -N 1

# specify number of tasks/cores per node required

# specify the walltime e.g 20 mins
#SBATCH -t 24:00:00

# set to email at start,end and failed jobs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robert.clampett@ucdconnect.ie

# run from current directory
cd $SLURM_SUBMIT_DIR
# command to use
hostname
for((c=1;c<11;c++))
do
	./SimpleSquare 20 800 3000 $c
	./SimpleSquare 40 800 3000 $c
	./SimpleSquare 60 800 3000 $c
done