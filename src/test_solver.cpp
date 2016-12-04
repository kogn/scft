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
     double alpha = configSettings.Read<double>("alpha");
     double beta = configSettings.Read<double>("beta");
     double kappa = configSettings.Read<double>("kappa");
     double tau = configSettings.Read<double>("tau");
     double chiN = configSettings.Read<double>("chiN");
     double domain[DIM];
     domain[0] = configSettings.Read<double>("domain0");
     domain[1] = configSettings.Read<double>("domain1");
     int n = configSettings.Read<int>("Steps_on_chain");
     int bw = configSettings.Read<int>("Band_width");
     int m[DIM];
     m[0] = configSettings.Read<int>("Grid_Size_x");
     m[1] = configSettings.Read<int>("Grid_Size_y");

     std::string output_filename = configSettings.Read<std::string>("Output_filename");
     std::string input_filename = configSettings.Read<std::string>("Input_filename");


  Iterator test(n,bw,m,alpha,beta,kappa,tau,chiN,domain);
  double * field = (double *)malloc(sizeof(double)*test.md*2);
  double * mu= (double *)malloc(sizeof(double)*test.md*2);

  int md = test.md;
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
    for(int i = 0; i<m[0]; i++)
      for(int j = 0; j<m[1]; j++)
      {
        mu[i*m[1]+j] = -2*cos(2*M_PI*i/m[0]) + 0.2*cos(2*M_PI*j/m[1]);
        mu[i*m[1]+j+md] = -mu[i*m[1]+j];
        field[i*m[1]+j] = mu[i*m[1]+j] - mu[i*m[1]+j+md];
        field[i*m[1]+j+md] = mu[i*m[1]+j] + mu[i*m[1]+j+md];
      }
  }

  test.solve_eqn(field);

  for(int i = 0; i<= test.n_step; i++)
    std::cout<<test.ptnfn(i)<<std::endl;
  free(field);
  free(mu);

  return 0;
}

