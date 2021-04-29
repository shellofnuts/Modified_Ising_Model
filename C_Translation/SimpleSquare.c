#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* Global Variable Declaration */
int anisotropyExchange = 0;
int anisotropyAxis = 1;

int initLatticeSize = 5;

int initWarmupSteps = 10;
int maxSimSteps = 100;
int minSimSteps = 50;

float startTemp = 0.4;
float endTemp = 1.6;
float intervalTemp = 1.05;
float criticalEstimate = 1.2;

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

int randomVec(double *output_vec){
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
    //printf("<%f, %f, %f> \n", output_vec[0], output_vec[1], output_vec[2]);
}

double getRandom(){
    return ((double)rand() / (double)RAND_MAX);
}

double gaussianRand(){
    double x,y,rsq,f;
    do {
        x = 2.0 * rand() / (double)RAND_MAX - 1.0;
        y = 2.0 * rand() / (double)RAND_MAX - 1.0;
        rsq = x * x + y * y;
    }while( rsq == 0.0 || rsq > 1.0 );

    f = sqrt( -2.0 * log(rsq) / rsq );

    return x*f;
}

double Hamiltonian(double *p0, double neighbours[nearestNeighbours][3]){
    double energy = 0;

    // Sum the altered dot product
    for(int i=0; i<nearestNeighbours;i++){
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

    double deltaE = originalE - newE;

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
    for(int i = 0; i < initLatticeSize; i++){
        for(int j = 0; j < initLatticeSize; j++){
            double neighbours[nearestNeighbours][3];
            double *p0, *p0_new;

            memset(neighbours, 0, nearestNeighbours*3*sizeof(double));
            p0 = malloc(3*sizeof(double));
            p0_new = malloc(3*sizeof(double));

            randomVec(p0_new);

            for(int k = 0; k < 3; k++){
                p0[k] = lattice[i][j][k];
                neighbours[0][k] = (i==0) ? lattice[initLatticeSize-1][j][k] : lattice[i-1][j][k];
                neighbours[1][k] = lattice[(i+1) % initLatticeSize][j][k];
                neighbours[2][k] = (j==0) ? lattice[i][initLatticeSize - 1][k] : lattice[i][j-1][k];
                neighbours[3][k] = lattice[i][(j+1) % initLatticeSize][k];
                }

            if(acceptChange(p0, p0_new, neighbours, T) == 1){
                printf("Changed: %d\n", j);
                for(int k = 0;k < 3; k++){
                    lattice[i][j][k] = p0_new[k];
                }
            }

            free(neighbours);
            free(p0);
            free(p0_new);
        }
    }
}

int warmup(double lattice[initLatticeSize][initLatticeSize][3], int maxSteps, double T){
    for(int i = 0; i< maxSteps; i++){
        alterLattice(lattice, T);
    }
}

int Metropolis(double lattice[initLatticeSize][initLatticeSize][3], int testSteps, double T, double *magSamples){
    for(int i = 0; i < testSteps; i++){
        alterLattice(lattice, T);
        magSamples[i] = magneticMoment(lattice);
    }
}

double magneticMoment(double lattice[initLatticeSize][initLatticeSize][3]){
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
    //printf("%f, %f, %f \n", avg_magVec[0], avg_magVec[1], avg_magVec[2]);
    return sqrt((avg_magVec[0]*avg_magVec[0])+(avg_magVec[1]*avg_magVec[1])+(avg_magVec[2]*avg_magVec[2]));
}

int main(void){
    int i, j, k;
    double t;

    FILE *fp;
    fp = fopen("10x10_spinDist.csv", "w+");

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

    for(i=0; i < initLatticeSize; i++){
        for(j=0; j < initLatticeSize; j++){
            printf("%f \t", lattice[i][j][2]);
        }
        printf("\n");
    }

    printf("%f\n", magneticMoment(lattice));
    warmup(lattice, initWarmupSteps, startTemp);
    for(i=0; i < initLatticeSize; i++){
        for(j=0; j < initLatticeSize; j++){
            printf("%f \t", lattice[i][j][2]);
        }
        printf("\n");
    }
    printf("%f\n", magneticMoment(lattice));
    /*t = startTemp;

    int simSteps = minSimSteps;

    while(t < endTemp){
            double sampleMags[simSteps];
            Metropolis(lattice, simSteps, t, sampleMags);
            fprintf(fp, "%f", t);
            for(i = 0; i < simSteps; i++){
                fprintf(fp, ", %f", sampleMags[i]);
            }
            fprintf(fp, "\n");
            t *= intervalTemp;
            free(sampleMags);
    }
    */
    free(lattice);
    fclose(fp);

}
