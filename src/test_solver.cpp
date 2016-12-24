#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "solver.h"
#include "homogeneous.hpp"
#include "Config.h"
#include <cmath>

int main(int argc, char * argv[])
{
    if(argc < 2)
    {
        std::cout<<"USE: ./main <path to config file>."<<std::endl;
        return 0;
    }
    Config configSettings(argv[1]);

    Homogeneous<Solver> test(configSettings);

    double * field = (double *)malloc(sizeof(double)*test.A.md*2);
    double * mu= (double *)malloc(sizeof(double)*test.A.md*2);

    int md = test.A.md;
    if(DIM == 1)
    {
      for(int i = 0; i<md; i++)
      {
        mu[i] = -2*cos(2*M_PI*i/md);
        mu[i+md] = -mu[i];
        field[i] = mu[i] - mu[i+md];
        field[i+md] = mu[i] + mu[i+md];
      }
    }
    if(DIM == 2)
    {
      for(int i = 0; i<test.A.m[0]; i++)
        for(int j = 0; j<test.A.m[1]; j++)
        {
          mu[i*test.A.m[1]+j] = -2*cos(2*M_PI*i/test.A.m[0]) + 0.2*cos(2*M_PI*j/test.A.m[1]);
          mu[i*test.A.m[1]+j+md] = -mu[i*test.A.m[1]+j];
          field[i*test.A.m[1]+j] = mu[i*test.A.m[1]+j] - mu[i*test.A.m[1]+j+md];
          field[i*test.A.m[1]+j+md] = mu[i*test.A.m[1]+j] + mu[i*test.A.m[1]+j+md];
        }
    }
    if(DIM == 3)
    {
        int * m = test.A.m;
        for(int i = 0; i<m[0]; i++)
            for(int j = 0; j<m[1]; j++)
                for(int k = 0; k<m[2]; k++)
                {
                    mu[i*m[1]*m[2]+j*m[2]+k] = -.2*cos(2*M_PI*i/m[0]) + 2*cos(2*M_PI*j/m[1]) + .2*cos(2*M_PI*k/m[2]);
                    mu[i*m[1]*m[2]+j*m[2]+k+md] = -mu[i*m[1]*m[2]+j*m[2]+k];
                    field[i*m[1]*m[2]+j*m[2]+k] = mu[i*m[1]*m[2]+j*m[2]+k] - mu[i*m[1]*m[2]+j*m[2]+k+md];
                    field[i*m[1]*m[2]+j*m[2]+k+md] = mu[i*m[1]*m[2]+j*m[2]+k] + mu[i*m[1]*m[2]+j*m[2]+k+md];
                }
    }

    test.density(field);


    for(int i = 0; i<= test.n_step; i++)
        std::cout<<test.ptnfn(i)<<std::endl;
    free(field);
    free(mu);

    return 0;
}

