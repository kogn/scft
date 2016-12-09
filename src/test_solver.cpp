#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "iterator.hpp"
#include "Config.h"

int main(int argc, char * argv[])
{
    if(argc < 2)
    {
        std::cout<<"USE: ./main <path to config file>."<<std::endl;
        return 0;
    }
    Config configSettings(argv[1]);

    Solver test(configSettings);

    double * field = (double *)malloc(sizeof(double)*test.md*2);
    double * mu= (double *)malloc(sizeof(double)*test.md*2);

    int md = test.md;
    /* if(DIM == 1) */
    /* { */
    /*   for(int i = 0; i<md; i++) */
    /*   { */
    /*     mu[i] = -2*cos(2*M_PI*i/md); */
    /*     mu[i+md] = -mu[i]; */
    /*     field[i] = mu[i] - mu[i+md]; */
    /*     field[i+md] = mu[i] + mu[i+md]; */
    /*   } */
    /* } */
    /* if(DIM == 2) */
    /* { */
    /*   for(int i = 0; i<m[0]; i++) */
    /*     for(int j = 0; j<m[1]; j++) */
    /*     { */
    /*       mu[i*m[1]+j] = -2*cos(2*M_PI*i/m[0]) + 0.2*cos(2*M_PI*j/m[1]); */
    /*       mu[i*m[1]+j+md] = -mu[i*m[1]+j]; */
    /*       field[i*m[1]+j] = mu[i*m[1]+j] - mu[i*m[1]+j+md]; */
    /*       field[i*m[1]+j+md] = mu[i*m[1]+j] + mu[i*m[1]+j+md]; */
    /*     } */
    /* } */
    test.density(field);


    for(int i = 0; i<= test.n_step; i++)
        std::cout<<test.ptnfn(i)<<std::endl;
    free(field);
    free(mu);

    return 0;
}

