#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "solver.h"
#include <cmath>

int main(int argc, char * argv[])
{
  double sum = 0;

  Solver test_solver(100,16,64,1,1,0,0,1);

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

  for(int j = 0; j<test_solver.md; j++)
    for(int i = 0; i<test_solver.n_coeff; i++)
    {
      sum += fabs(test_solver.spectdata[test_solver.n_coeff*j+i][0] - j);
      sum += fabs(test_solver.spectdata[test_solver.n_coeff*j+i][1] - i);
    }
  std::cout<<sum<<std::endl;

  return 0;
}

