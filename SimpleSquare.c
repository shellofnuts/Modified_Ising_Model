#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/* Defined Functions */

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

/* Global Variable Declaration */
double anisotropyExchange = 0.;
double anisotropyAxis = 0.;

int initLatticeSize;

int initWarmupSteps = 200;
int maxSimSteps = 1000;
int minSimSteps = 300;

double criticalEstimate = 1.2;
double startTemp;
double endTemp;
double intervalTemp = 1.0005;

int simNumber = 0;
int nearestNeighbours = 4;

/* Function Declaration */
double gaussianRand();
int randomVec(double *output_vec);
double getRandom();
double Hamiltonian(double *p0, double neighbours[nearestNeighbours][3]);
double acceptanceRatio(double energy, double T);
int acceptChange(double *p0, double *p0_new, double neighbours[nearestNeighbours][3], double Tk);
int alterLattice(double lattice[initLatticeSize][initLatticeSize][3], double T);
int warmup(double lattice[initLatticeSize][initLatticeSize][3], int maxSteps, double T);
double magneticMoment(double lattice[initLatticeSize][initLatticeSize][3]);
int Metropolis(double lattice[initLatticeSize][initLatticeSize][3], int testSteps, double T, double *magSamples);
int updateSimStep(double T, double scale);

/* Function Description */

int updateSimStep(double T, double scale){
    /* Inputs T.
       Returns the number of simsteps associated with that T.
       This is calculated by scaling the simsteps to a x^2 centred on criticalEstimate */

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
        output_vec[i] = new_vec[i]/vec_norm;
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

double Hamiltonian(double *p0, double neighbours[nearestNeighbours][3]){
    /*
        Inputs pointers to point of interest and an array of NearestNeighbours
        Calculates and returns energy associated.
    */

    double energy = 0;
	int i;

    // Sum the altered dot product
    for(i=0; i<nearestNeighbours;i++){
        energy += -1*((1-anisotropyExchange)*(p0[0]*neighbours[i][0] + p0[1]*neighbours[i][1]) + p0[2]*neighbours[i][2]);
    }
    // Add the easy axis anisotropy
    energy += -2*anisotropyAxis*(p0[2]*p0[2]);
    return energy;
}

double acceptanceRatio(double energy, double T){
    return exp(-1*energy/T);
}

int acceptChange(double *p0, double *p0_new, double neighbours[nearestNeighbours][3], double Tk){
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

int alterLattice(double lattice[initLatticeSize][initLatticeSize][3], double T){
	int i, j, k;
    for(i = 0; i < initLatticeSize; i++){
        for(j = 0; j < initLatticeSize; j++){
            double neighbours[nearestNeighbours][3];
            double *p0, *p0_new;

            memset(neighbours, 0, nearestNeighbours*3*sizeof(double));
            p0 = malloc(3*sizeof(double));
            p0_new = malloc(3*sizeof(double));

            randomVec(p0_new);

            for(k = 0; k < 3; k++){
                p0[k] = lattice[i][j][k];
                neighbours[0][k] = (i==0) ? lattice[initLatticeSize-1][j][k] : lattice[i-1][j][k];
                neighbours[1][k] = lattice[(i+1) % initLatticeSize][j][k];
                neighbours[2][k] = (j==0) ? lattice[i][initLatticeSize - 1][k] : lattice[i][j-1][k];
                neighbours[3][k] = lattice[i][(j+1) % initLatticeSize][k];
                }

            if(acceptChange(p0, p0_new, neighbours, T) == 1){
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

int warmup(double lattice[initLatticeSize][initLatticeSize][3], int maxSteps, double T){
	int i;
    for(i = 0; i< maxSteps; i++){
        alterLattice(lattice, T);
    }
}

int Metropolis(double lattice[initLatticeSize][initLatticeSize][3], int testSteps, double T, double *magSamples){
    /*
        Loops through a set of simulation steps, altering the lattice when energy allows.
        It saves the average magnetisation of the lattice to output pointer array magSamples.
    */

	int i;

    for(i = 0; i < testSteps; i++){
        alterLattice(lattice, T);
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
    int i, j, k;
    double t;

	opterr = 0;
	int c;

	char options[] = "N:U:L:i:a:e:T:";

	while((c = getopt(argc, argv, options)) != -1){
        switch (c)
        {
        case 'N':
            initLatticeSize = atoi(optarg);
            break;
        case 'U':
            maxSimSteps = atoi(optarg);
            break;
        case 'L':
            minSimSteps = atoi(optarg);
            break;
        case 'i':
            simNumber = atoi(optarg);
            break;
        case 'a':
            anisotropyAxis = atof(optarg);
            break;
        case 'e':
            anisotropyExchange = atof(optarg);
            break;
        case 'T':
            criticalEstimate = atof(optarg);
            break;
        case '?':
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
    // Print to console the initialised variables for future reference.
    printf("Initialised with:\n");
    printf("Lattice Size = %d\n", initLatticeSize);
    printf("Maximum Simulation Steps = %d\n", maxSimSteps);
    printf("Minimum Simulation Steps = %d\n", minSimSteps);
    printf("Axis Anisotropy Value: %f\n", anisotropyAxis);
    printf("Exchange Anisotropy Value: %f\n", anisotropyExchange);
    printf("Critical Estimate: %f\n", criticalEstimate);

    FILE *fp;

    char filename[32];
    sprintf(filename, "%dx%d_spinDist_%d.csv", initLatticeSize, initLatticeSize, simNumber);

    fp = fopen(filename, "w+");

	// Update Global Variables
	startTemp = criticalEstimate - 0.15;
	endTemp = criticalEstimate + 0.15;

	// Set a random seed for random numbers
	srand(time(NULL));

    // Initialise the Lattice

    double lattice[initLatticeSize][initLatticeSize][3];
    memset(lattice, 0, initLatticeSize*initLatticeSize*sizeof(double));

    for(i=0; i<initLatticeSize; i++){
        for(j=0; j<initLatticeSize; j++){
                lattice[i][j][0] = 0;
                lattice[i][j][1] = 0;
                lattice[i][j][2] = 1;
        }
    }

    warmup(lattice, initWarmupSteps, startTemp);

    t = startTemp;

    int simSteps = minSimSteps;
    double scale = max((criticalEstimate - startTemp), (endTemp - criticalEstimate));

    while(t < endTemp){
            double sampleMags[simSteps];
            memset(sampleMags, 0, simSteps*sizeof(double));
            Metropolis(lattice, simSteps, t, sampleMags);
            fprintf(fp, "%f", t);
            for(i = 0; i < simSteps; i++){
                fprintf(fp, ", %f", sampleMags[i]);
            }
            fprintf(fp, "\n");
            t *= intervalTemp;
            simSteps = updateSimStep(t, scale);
            //free(sampleMags);
    }
    //free(lattice);
    fclose(fp);

	return 0;
}
