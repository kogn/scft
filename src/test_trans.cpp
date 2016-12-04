#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "solver.h"
#include <cmath>
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


  double sum = 0;

  Solver test_solver(n,bw,m,alpha,beta,kappa,tau,domain);

  for(int j = 0; j<test_solver.md; j++)
    for(int i = 0; i<test_solver.n_coeff; i++)
    {
      test_solver.spectdata[test_solver.n_coeff*j+i][0] = j;
      test_solver.spectdata[test_solver.n_coeff*j+i][1] = i;
    }
  test_solver.inv_so3();
  test_solver.inv_space();
  test_solver.for_space();
  test_solver.for_so3();

  test_solver.inv_so3();
  test_solver.inv_space();
  test_solver.for_space();
  test_solver.for_so3();

  test_solver.inv_so3();
  test_solver.inv_space();
  test_solver.for_space();
  test_solver.for_so3();

  for(int j = 0; j<test_solver.md; j++)
    for(int i = 0; i<test_solver.n_coeff; i++)
    {
      sum += fabs(test_solver.spectdata[test_solver.n_coeff*j+i][0] - j);
      sum += fabs(test_solver.spectdata[test_solver.n_coeff*j+i][1] - i);
    }
  std::cout<<sum<<std::endl;

  return 0;
}

