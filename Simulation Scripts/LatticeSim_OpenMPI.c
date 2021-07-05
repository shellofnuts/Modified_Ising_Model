#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "mpi.h"
#include <stddef.h>

/* Defined Macros */

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

/* Structure Definition */

typedef struct {
    int apos;
    int bpos;
    int FM_type;
}   NNvec_t;

/*
*   Struct stores information on direction of nearest neighbour
*   and the type of bond interaction (FM: +1, AFM: -1).
*   This is distributed in MPI_Bcast in main().
*
*   Used in readFile(), alterLattice(), Warmup() & Metropolis().
*/

/* Global Variable Declaration */

double anisotropyExchange = 0.;         // Default values for magnetic anisotropy.
double anisotropyAxis = 0.;             // Values are altered by command line arguments.

int initLatticeSize;                    // No default value, must specify in CMD line args.

int initWarmupSteps = 500;
int maxSimSteps = 1000;
int minSimSteps = 300;

double criticalEstimate = 101.0;        // Estimate for Tc. Centres the temperature range search.

double startTemp;                       // Local variables for each MPI node.
double endTemp;

double intervalTemp = 5.0;
double start;

double rangeTemp = 200.0;

int simNumber = 0;
int nearestNeighbours = 4;

/* Constants */

const double boltzmann = 8.6173E-5; 	// Boltzmann Constant in eV
const double J_bond = 0.86E-3; 			// J interaction energy in eV
const double uB_moment = 3.0;			// Intrinsic Lattice Moment

/* Function Declaration */

double gaussianRand();
int randomVec(double *output_vec);
double getRandom();
double Hamiltonian( double *p0, double neighbours[nearestNeighbours][3], int FM_type[nearestNeighbours]);
double acceptanceRatio(double energy, double T);
int acceptChange(double *p0, double *p0_new, double neighbours[nearestNeighbours][3], int FM_type[nearestNeighbours], double Tk);
int alterLattice(double lattice[initLatticeSize][initLatticeSize][3], double T, NNvec_t set_of_NN[nearestNeighbours]);
int warmup(double lattice[initLatticeSize][initLatticeSize][3], int maxSteps, double T, NNvec_t set_of_NN[nearestNeighbours]);
double magneticMoment(double lattice[initLatticeSize][initLatticeSize][3]);
int Metropolis(double lattice[initLatticeSize][initLatticeSize][3], int testSteps, double T, double *magSamples, NNvec_t set_of_NN[nearestNeighbours]);
int updateSimStep(double T, double scale);
int readFile(FILE *input_file, NNvec_t vectors[nearestNeighbours]);

/* Function Description */

int readFile(FILE *input_file, NNvec_t vectors[nearestNeighbours]){
    // Reads INPUTVEC file.
    // Extracts the nearest neighbour directions and the magnetic
    // interaction type of each NN interaction.

    size_t count = 0;
    while(fscanf(input_file, "%d,%d,&d", &vectors[count].apos, &vectors[count].bpos, &vectors[count].FM_type) == 3){
        count++;
    }
}

int updateSimStep(double T, double scale){
    /* Inputs T.
       Returns the number of simsteps associated with that T.
       This is calculated by scaling the simsteps to a x^2 centred on criticalEstimate.
    */
    double steps = maxSimSteps - (maxSimSteps-minSimSteps)*((T-criticalEstimate)*(T-criticalEstimate)/(scale*scale));
    int result = (int) steps;
    //printf("%f, %d\n", steps, result);
    return result;
}

int randomVec(double *output_vec){
    /*
        Input is a pointer to a 1D vector array.
        Creates 3 Gaussian random numbers, then normalises.
        Assignes to output pointer.
    */

    int i;
    double new_vec[3] = {0};
    for(i = 0; i<3;i++){
        new_vec[i] = gaussianRand();
    }
    double vec_sum = 0;
    for(i=0;i<3;i++){
        vec_sum += new_vec[i]*new_vec[i];
    }
    double vec_norm = sqrt(vec_sum);
    for(i=0;i<3;i++){
        output_vec[i] = uB_moment * new_vec[i] / vec_norm;
    }
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
    double x,y,rsq,f;
    do {
        x = 2.0 * (rand() / (double)RAND_MAX) - 1.0;
        y = 2.0 * rand() / (double)RAND_MAX - 1.0;
        rsq = x * x + y * y;
    }while( rsq == 0.0 || rsq > 1.0 );
    f = sqrt( -2.0 * log(rsq) / rsq );
    return x*f;
}

double Hamiltonian( double *p0, double neighbours[nearestNeighbours][3], int FM_type[nearestNeighbours]){
    /*
        Inputs pointers to spin, array of NearestNeighbours and interaction type.
        Calculates and returns energy associated.
    */

    double energy = 0.0;
	int i;

    // Sum the altered dot product
    for(i=0; i<nearestNeighbours;i++){
        energy += FM_type[i] * ((1 - anisotropyExchange)*(p0[0]*neighbours[i][0] + p0[1]*neighbours[i][1]) + p0[2]*neighbours[i][2]);
    }

	energy *= -1 * J_bond;

    // Add the easy axis anisotropy
    energy += -1 * anisotropyAxis * ( p0[2] * p0[2] );
    return energy;
}

double acceptanceRatio(double energy, double T){
	double Tk = boltzmann * T;		// Convert T to "real" temperature instead of normalised T (T/kb).
    return exp(-1 * energy / Tk);
}

int acceptChange(double *p0, double *p0_new, double neighbours[nearestNeighbours][3], int FM_type[nearestNeighbours], double Tk){
    double originalE = Hamiltonian(p0, neighbours, FM_type);
    double newE = Hamiltonian(p0_new, neighbours, FM_type);

    double deltaE = newE - originalE;
    /*
    There are three cases:
    1. If the energy is decreased by change, accept change.
    2. If the energy increases but is acceptable due to temp, accept change.
    3. If the energy increases but is not acceptable due to temp.

    We only need to check the first two cases as they are the only ones that make a change.
    */

    if(deltaE <= 0){
        return 1;
    }
    else if(getRandom() <= acceptanceRatio(deltaE, Tk)){
        return 1;
    }
    else{
        return 2;
    }
}

int alterLattice(double lattice[initLatticeSize][initLatticeSize][3], double T, NNvec_t set_of_NN[nearestNeighbours]){
	int i, j, k, l;
    for(i = 0; i < initLatticeSize; i++){
        for(j = 0; j < initLatticeSize; j++){
            double neighbours[nearestNeighbours][3];
            double *p0, *p0_new;
            int *FM_type;

            memset(neighbours, 0, nearestNeighbours*3*sizeof(double));
            p0 = malloc(3*sizeof(double));
            p0_new = malloc(3*sizeof(double));
            FM_type = malloc(nearestNeighbours*sizeof(int));

            int A, B;

            randomVec(p0_new);

            for(k = 0; k < 3; k++){
                p0[k] = lattice[i][j][k];
            }
            for(l = 0; l< nearestNeighbours; l++){
                    /*
                     *   There are four scenarios that need to be dealt with for boundary conditions:
                     *
                     *  1. A >= 0 and B >= 0, this just needs modulo operator.
                     *  2. A = -1 and i = 0, then we need set the array position to the last position on i axis
                     *  3. B = -1 and j = 0, same issue.
                     *  4. A and B = -1 and i and j = 0, combination of 2 and 3.
                     *
                     *  There may be a way to combine scenarios 2-4 but currently not sure how to do that.
                     *  Ideally, I would only do a quick check to see if i or j are 0 first, so that I do less if-else computations.

                     *  Unpack the
                     */
                A = (i == 0 && set_of_NN[l].apos < 0) ? (initLatticeSize -1) : (i + set_of_NN[l].apos) % initLatticeSize;
                B = (j == 0 && set_of_NN[l].bpos < 0) ? (initLatticeSize -1) : (i + set_of_NN[l].bpos) % initLatticeSize;
                for(k = 0; k < 3; k++){
                    neighbours[l][k] = lattice[A][B][k];
                }

                FM_type[l] = set_of_NN[l].FM_type;

            }

            if(acceptChange(p0, p0_new, neighbours, FM_type, T) == 1){
                for(k = 0;k < 3; k++){
                    lattice[i][j][k] = p0_new[k];
                }
            }

            //free(neighbours);
            free(p0);
            free(p0_new);
        }
    }
}

int warmup(double lattice[initLatticeSize][initLatticeSize][3], int maxSteps, double T, NNvec_t set_of_NN[nearestNeighbours]){
	int i;
    for(i = 0; i< maxSteps; i++){
        alterLattice(lattice, T, set_of_NN);
    }
}

int Metropolis(double lattice[initLatticeSize][initLatticeSize][3], int testSteps, double T, double *magSamples, NNvec_t set_of_NN[nearestNeighbours]){
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

double magneticMoment(double lattice[initLatticeSize][initLatticeSize][3]){
    /*
        Inputs lattice array.
        Calculates the average magnetic moment.
        Returns magnitude of vector.
    */

    int i, j;
    double avg_magVec[3] = {0};
    for(i = 0; i < initLatticeSize; i++){
        for(j = 0; j < initLatticeSize; j++){
            avg_magVec[0] += lattice[i][j][0];
            avg_magVec[1] += lattice[i][j][1];
            avg_magVec[2] += lattice[i][j][2];
        }
    }
    for(i =0; i < 3; i++){
        avg_magVec[i] /= (double)(initLatticeSize*initLatticeSize);
    }
    return sqrt((avg_magVec[0]*avg_magVec[0])+(avg_magVec[1]*avg_magVec[1])+(avg_magVec[2]*avg_magVec[2]));
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
    /*
        Implement command line arguments using getopt().
        char options[] stores the command line letter. ":" denotes that a value after the letter is required.
    */
	char options[] = "N:U:L:i:a:e:T:s:r:";

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

        case 'a':
            // Set axis anisotropy value
            anisotropyAxis = atof(optarg);
            break;

        case 'e':
            // Set exchange anisotropy value
            anisotropyExchange = atof(optarg);
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

		case 'r':
			// Set the Temperature range of the system.
			rangeTemp = atof(optarg);
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
        printf("Axis Anisotropy Value: %f\n", anisotropyAxis);
        printf("Exchange Anisotropy Value: %f\n", anisotropyExchange);
		printf("Magnetic Interaction Energy: %f\n \n", J_bond);
        printf("Critical Estimate: %f\n", criticalEstimate);
		printf("Temperature Range: %f\n", rangeTemp);
    }

	// Set a random seed for random numbers
	srand(time(NULL));

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
        printf("Distributing NNpos...\n");
    }

    /*
     *  Set up the MPI_Datatype struct to be able to distribute the struct data to all nodes.
     *  We need to define the number of elements per message, type of elements and their memory offsets.
     *  We create a temporary type and then calculate the excess memory used in each memory.
     *  Output type is then a resized tmp type.
    */

	int elements = 3;
	int array_of_blocklengths[] = {1, 1, 1};
	MPI_Datatype array_of_types[] = {MPI_INT, MPI_INT, MPI_INT};
	MPI_Aint array_of_displacements[] = { offsetof(NNvec_t, apos), offsetof(NNvec_t, bpos), offsetof(NNvec_t, FM_type)};
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

    // Set the output name here for each individual process.
    // Will use process 0 to collate files into one file.

    sprintf(filename, "%dx%d_spinDist_%d_%d.csv", initLatticeSize, initLatticeSize, simNumber, myrank);

    fp = fopen(filename, "w+");

	// Update Global Variables according to process.
	double span = rangeTemp / commsize;
	startTemp = criticalEstimate - rangeTemp/2 + myrank*span;
	endTemp = startTemp + span;

    // Initialise the Lattice
    double lattice[initLatticeSize][initLatticeSize][3];
    memset(lattice, 0, initLatticeSize*initLatticeSize*sizeof(double));

    for(i=0; i<initLatticeSize; i++){
        for(j=0; j<initLatticeSize; j++){
                lattice[i][j][0] = 0.0;
                lattice[i][j][1] = 0.0;
                lattice[i][j][2] = uB_moment;
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
            t += intervalTemp;
            //simSteps = updateSimStep(t, scale);
    }
    fclose(fp);

    MPI_Barrier(MPI_COMM_WORLD);

    if(myrank==0){
        printf("Execution Time: %f\n", MPI_Wtime() - start);
    }
    MPI_Finalize();

	return 0;
}
