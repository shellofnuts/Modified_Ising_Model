#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "mpi.h"
#include <stddef.h>

/* Defined Functions */

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

/* Structure Definition */

typedef struct {
    int apos;
    int bpos;
}   NNvec_t;

/* Global Variable Declaration */

int initLatticeSize;

int initWarmupSteps = 500;
int maxSimSteps = 1000;
int minSimSteps = 300;

double criticalEstimate = 2.2;
double startTemp;
double endTemp;
double intervalTemp = 1.0005;
double start;

int states[4] = {3, 1, -1, -3};

int simNumber = 0;
int nearestNeighbours = 4;

/* Constants */

const double boltzmann = 8.6173E-28; // Boltzmann Constant in eV
const double J_bond = 0.86E-3; // J interaction energy in eV

/* Function Declaration */
double gaussianRand();
int randomVec();
double getRandom();
int Hamiltonian(int p0, int neighbours[nearestNeighbours]);
double acceptanceRatio(int energy, double T);
int acceptChange(int p0, int p0_new, int neighbours[nearestNeighbours], double Tk);
int alterLattice(int lattice[initLatticeSize][initLatticeSize], double T, NNvec_t set_of_NN[nearestNeighbours]);
int warmup(int lattice[initLatticeSize][initLatticeSize], int maxSteps, double T, NNvec_t set_of_NN[nearestNeighbours]);
double magneticMoment(int lattice[initLatticeSize][initLatticeSize]);
int Metropolis(int lattice[initLatticeSize][initLatticeSize], int testSteps, double T, double *magSamples, NNvec_t set_of_NN[nearestNeighbours]);
int updateSimStep(double T, double scale);
int readFile(FILE *input_file, NNvec_t vectors[nearestNeighbours]);

/* Function Description */

int readFile(FILE *input_file, NNvec_t vectors[nearestNeighbours]){
    // Read NN directions
    size_t count = 0;
    while(fscanf(input_file, "%d,%d", &vectors[count].apos, &vectors[count].bpos) == 2){
        count++;
    }
}

int updateSimStep(double T, double scale){
    /* Inputs T.
       Returns the number of simsteps associated with that T.
       This is calculated by scaling the simsteps to a x^2 centred on criticalEstimate */

    double steps = maxSimSteps - (maxSimSteps-minSimSteps)*((T-criticalEstimate)*(T-criticalEstimate)/(scale*scale));
    int result = (int) steps;
    //printf("%f, %d\n", steps, result);
    return result;
}

int randomVec(){
    /*
        Returns a random value from the states set.
    */

    int i;
    int new_vec = gaussianRand();

    return new_vec;
}

double getRandom(){
    return rand() / (double)RAND_MAX;
}

double gaussianRand(){
    /*
        Uses the Box-Muller transform on a random uniform [0,1] to
        create a random gaussian distribution.
        Returns 1 number.
    */

    double x = 2 * getRandom();
    int output;

    if(x <= 1){
        output = states[0]; // +3
    }
    else if(x <= 2){
        output = states[1]; // +1
    }
    else if(x <= 3){
        output = states[2]; // -1
    }
    else{
        output = states[3]; // -3
    }
    return output;
}

int Hamiltonian(int p0, int neighbours[nearestNeighbours]){
    /*
        Inputs pointers to point of interest and an array of NearestNeighbours
        Calculates and returns energy associated.
    */

    double energy = 0;
	int i;

    // Sum the altered dot product
    for(i=0; i<nearestNeighbours;i++){
        energy += (p0*neighbours[i]);
    }
    return -1 * J_bond * energy;
}

double acceptanceRatio(int energy, double T){
    double Tk = boltzmann * T;
    return exp(-1*energy/Tk);
}

int acceptChange(int p0, int p0_new, int neighbours[nearestNeighbours], double Tk){
    double originalE = Hamiltonian(p0, neighbours);
    double newE = Hamiltonian(p0_new, neighbours);

    double deltaE = newE - originalE;
    /*
    There are three cases:
    1. If the energy is decreased by change, accept change.
    2. If the energy increases but is acceptable due to temp, accept change.
    3. If the energy increases but is not acceptable due to temp.

    We only need to check the first two cases as they are the only ones that make a change.
    */

    if(deltaE < 0){
        return 1;
    }
    else if(getRandom() < acceptanceRatio(deltaE, Tk)){
        return 1;
    }
    else{
        return 2;
    }
}

int alterLattice(int lattice[initLatticeSize][initLatticeSize], double T, NNvec_t set_of_NN[nearestNeighbours]){
	int i, j, k, l;
    for(i = 0; i < initLatticeSize; i++){
        for(j = 0; j < initLatticeSize; j++){
            int neighbours[nearestNeighbours];
            int p0, p0_new;

            memset(neighbours, 0, nearestNeighbours*sizeof(int));

            int A, B;

            p0_new = randomVec();
            p0 = lattice[i][j];

            for(l = 0; l< nearestNeighbours; l++){
                    /*
                        There are four scenarios that need to be dealt with for boundary conditions:

                        1. A >= 0 and B >= 0, this just needs modulo operator.
                        2. A = -1 and i = 0, then we need set the array position to the last position on i axis
                        3. B = -1 and j = 0, same issue.
                        4. A and B = -1 and i and j = 0, combination of 2 and 3.

                        There may be a way to combine scenarios 2-4 but currently not sure how to do that.
                        Ideally, I would only do a quick check to see if i or j are 0 first, so that I do less if-else computations.
                    */
                    A = (i == 0 && set_of_NN[l].apos < 0) ? (initLatticeSize -1) : (i + set_of_NN[l].apos) % initLatticeSize;
                    B = (j == 0 && set_of_NN[l].bpos < 0) ? (initLatticeSize -1) : (i + set_of_NN[l].bpos) % initLatticeSize;
                    neighbours[l] = lattice[A][B];
            }

            if(acceptChange(p0, p0_new, neighbours, T) == 1){
                lattice[i][j] = p0_new;
            }
        }
    }
}

int warmup(int lattice[initLatticeSize][initLatticeSize], int maxSteps, double T, NNvec_t set_of_NN[nearestNeighbours]){
	int i;
    for(i = 0; i< maxSteps; i++){
        alterLattice(lattice, T, set_of_NN);
    }
}

int Metropolis(int lattice[initLatticeSize][initLatticeSize], int testSteps, double T, double *magSamples, NNvec_t set_of_NN[nearestNeighbours]){
    /*
        Loops through a set of simulation steps, altering the lattice when energy allows.
        It saves the average magnetisation of the lattice to output pointer array magSamples.
    */

	int i;

    for(i = 0; i < testSteps; i++){
        alterLattice(lattice, T, set_of_NN);
        magSamples[i] = magneticMoment(lattice);
    }
}

double magneticMoment(int lattice[initLatticeSize][initLatticeSize]){
    /*
        Inputs lattice array.
        Calculates the average magnetic moment.
        Returns magnitude of vector.
    */

    int i, j;
    double avg_magVec = 0;
    for(i = 0; i < initLatticeSize; i++){
        for(j = 0; j < initLatticeSize; j++){
            avg_magVec += lattice[i][j];
        }
    }
    avg_magVec /= (double)(initLatticeSize*initLatticeSize);
    return avg_magVec;
}

int main(int argc, char *argv[]){
    /*
    The approach taken here is to multithread the temperature average section.
    Each process will take a certain chunk of the temperature range and perform the calculations necessary.
    This parallisation of the temperature sweep function should hopefully decrease sim times.
    */

    int myn, myrank, commsize;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize); MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int i, j, k;
    double t;

	opterr = 0;
	int c;

	char options[] = "N:U:L:i:a:e:T:s:";

	while((c = getopt(argc, argv, options)) != -1){
        switch (c)
        {
        case 'N':
            // Set Lattice Size
            initLatticeSize = atoi(optarg);
            break;

        case 'U':
            // Set Upper Bound
            maxSimSteps = atoi(optarg);
            break;

        case 'L':
            // Set Lower Bound
            minSimSteps = atoi(optarg);
            break;

        case 'i':
            // Set iteration number
            simNumber = atoi(optarg);
            break;

        case 'T':
            // Set estimate of critcal temperature
            criticalEstimate = atof(optarg);
            break;

        case 's':
            // Set the symmetry of the system
            // Required for non-square lattices
            nearestNeighbours = atoi(optarg);
            break;

        case '?':
            // Exception handling
            if (strchr(options, optopt) != NULL){
                printf("Option -%c requires an argument.\n", optopt);
            }
            else if (isprint(optopt)){
                printf ("Unknown option `-%c'.\n", optopt);
            }
            else{
                printf("Unknown option character `\\x%x'.\n", optopt);
            }
            return 1;

        default:
            abort();
        }
	}
    if(myrank==0){
    // Print to console the initialised variables for future reference.
    printf("Initialised with:\n");
    printf("Lattice Size = %d\n", initLatticeSize);
    printf("Maximum Simulation Steps = %d\n", maxSimSteps);
    printf("Minimum Simulation Steps = %d\n", minSimSteps);
    printf("Critical Estimate: %f\n", criticalEstimate);
    }

    // Declare NN before reading file.
    NNvec_t nearestNeighbourPos[nearestNeighbours];

    /*
    Declare input file that holds symmetry vectors.
    Will be a set of 4 tuples for square, 6 for triangular etc.
    Currently need to declare nearestNeighbours on command line.
    Main process will read the file and then broadcast the data to all the rest of the process.
    */
    if(myrank==0){
        FILE *input_fp;
        input_fp = fopen("INPUTVECS", "r");
        readFile(input_fp, nearestNeighbourPos);
        fclose(input_fp);
    }
    /*
        Set up the MPI_Datatype struct to be able to distribute the struct data to all nodes.
    */

	int elements = 2;
	int array_of_blocklengths[] = {1, 1};
	MPI_Datatype array_of_types[] = {MPI_INT, MPI_INT};
	MPI_Aint array_of_displacements[] = { offsetof(NNvec_t, apos), offsetof(NNvec_t, bpos)};
	MPI_Datatype tmp_type, my_mpi_struct_type;
	MPI_Aint lb, extent;

	MPI_Type_create_struct(elements, array_of_blocklengths, array_of_displacements, array_of_types, &tmp_type);
	MPI_Type_get_extent(tmp_type, &lb, &extent);
	MPI_Type_create_resized(tmp_type, lb, extent, &my_mpi_struct_type);
	MPI_Type_commit(&my_mpi_struct_type);

	// Broadcast the data from process 0 to all other processes
    MPI_Bcast(nearestNeighbourPos, nearestNeighbours, my_mpi_struct_type, 0, MPI_COMM_WORLD);

    // Declare output csv file.
    FILE *fp;

    char filename[32];

    /*
        Set the output name here for each individual process.
        Will use process 0 to colate files into one file.
    */
    sprintf(filename, "%dx%d_spinDist_%d_%d.csv", initLatticeSize, initLatticeSize, simNumber, myrank);

    fp = fopen(filename, "w+");

	// Update Global Variables according to process.
	double span = 2.0 / commsize;
	startTemp = criticalEstimate - 1.0 + myrank*span;
	endTemp = startTemp + span;

	// Set a random seed for random numbers
	srand(time(NULL));

    // Initialise the Lattice
    int lattice[initLatticeSize][initLatticeSize];
    memset(lattice, 0, initLatticeSize*initLatticeSize*sizeof(int));

    for(i=0; i<initLatticeSize; i++){
        for(j=0; j<initLatticeSize; j++){
                lattice[i][j] = 3;
        }
    }

    // "Warmup" the lattice.
    warmup(lattice, initWarmupSteps, startTemp, nearestNeighbourPos);

    t = startTemp;

    // NOTE: Depreciating updateSimSteps func at the moment. Setting simsteps to Upper Bound (-U)

    int simSteps = maxSimSteps;
    //int simSteps = minSimSteps;
    //double scale = max((criticalEstimate - startTemp), (endTemp - criticalEstimate));

    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank==0){
        start = MPI_Wtime();
    }

    while(t < endTemp){
            double sampleMags[simSteps];
            memset(sampleMags, 0, simSteps*sizeof(double));
            Metropolis(lattice, simSteps, t, sampleMags, nearestNeighbourPos);
            fprintf(fp, "%f", t);
            for(i = 0; i < simSteps; i++){
                fprintf(fp, ", %f", sampleMags[i]);
            }
            fprintf(fp, "\n");
            t *= intervalTemp;
            //simSteps = updateSimStep(t, scale);
    }
    fclose(fp);

    MPI_Barrier(MPI_COMM_WORLD);

    if(myrank==0){
        printf("Time: %f\n", MPI_Wtime() - start);
        /*
        FILE *output_fp;
        char fileout[32];
        sprintf(fileout, "%dx%d_spinDist_%d.csv", initLatticeSize, initLatticeSize, simNumber);

        output_fp = fopen(filename, "w+");

        char buff[1024];
        for(i=0; i<commsize; i++){

        }

        */
    }
    MPI_Finalize();

	return 0;
}
