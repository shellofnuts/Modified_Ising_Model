#!/bin/bash -l
#SBATCH --job-name=CSQ_MAE
# speficity number of nodes
#SBATCH -n 8

# specify the walltime e.g 20 mins
#SBATCH -t 24:00:00

# run from current directory
cd $SLURM_SUBMIT_DIR

module load intel/intel-mpi/2020.4
module load intel/intel-ipp/2020.4
module load intel/intel-cc/2020.4
module load intel/intel-mkl/2020.4
module load openmpi/4.0.1


export PATH=/opt/software/openmpi/openmpi-4.0.1/bin:$PATH
export LD_LIBRARY_PATH=/opt/software/intel/2020u4/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin/

# command to use
# Write stdout+stderr to file
#SBATCH -o output.txt

ulimit -s unlimited

mpiexec hostname

rm ScriptLogs

for s in 40 ; do
for a in $(seq 0.7 0.1 1.0) ; do
	for((c=1; c<5; c++))
	do
		echo "Lattice Size: $s" >> ScriptLogs
		echo "MAE: $a" >> ScriptLogs
		mpirun -np 8 /your/path/here/to/script/LatticeSim_OpenMPI -N "$s" -U 100000 -i "$c" -a "$a" -T 12.7 -r 25.0 -I 0.5 > MyLogFile
		for((i=0; i<8; i++))
		do
			cat "$s"x"$s"_spinDist_"$c"_"$i".csv >> "$s"x"$s"_spinDistHex_a-"$a"_"$c".csv
			rm "$s"x"$s"_spinDist_"$c"_"$i".csv
		done
	done
done
done

module load anaconda

rm PythonLogFile

for a in $(seq 0.7 0.1 1.0) ; do
	echo "Parsing $a anisotropy data" >> PythonLogFile
	python Analysis.py -n 4 -o "HexOutput" -i "[20,40]" -a $a
done