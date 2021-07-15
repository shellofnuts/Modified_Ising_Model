#!/bin/bash -l
#SBATCH --job-name=CSQ_OpenMPI
# speficity number of nodes
#SBATCH -n 4

# specify the walltime e.g 20 mins
#SBATCH -t 10:00:00

# run from current directory
cd $SLURM_SUBMIT_DIR

module load intel/intel-mpi/2020.4
module load intel/intel-ipp/2020.4
module load intel/intel-cc/2020.4
module load intel/intel-mkl/2020.4
module load openmpi/4.0.1

which mpirun

#export PATH=/opt/software/openmpi/openmpi-4.0.1/bin:$PATH
export LD_LIBRARY_PATH=/opt/software/intel/2020u4/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin/

# command to use
# Write stdout+stderr to file
#SBATCH -o output.txt

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/software/intel/2019Parallel/lib/intel64_lin:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/software/intel/2019Parallel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64_lin:$LD_LIBRARY_PATH
ulimit -s unlimited

mpiexec hostname

mpirun -np 4 /home/people/20204902/SimpleSquare/LatticeSim_OpenMPI -N 20 -L 500 -U 1000 -a 0.8 -T 1.0 -i 1 > MyLogFile
for((i=1; i<5; i++))
do
	echo cat 20x20_spinDist_1_"$i".SEP > 20x20_spinDist_1.csv
done