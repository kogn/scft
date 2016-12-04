#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "iterator.hpp"

int main(int argc, char * argv[])
{
  Iterator test(200,8,64,10000000,10000000,5,5,6,2);
  double * field = (double *)malloc(sizeof(double)*test.md*2);
  double * mu= (double *)malloc(sizeof(double)*test.md*2);
  for(int i = 0;i<test.md; i++)
  {
    mu[i] = -2*cos(2*M_PI*i/test.md);
    mu[i+test.md] = -mu[i];
    field[i] = mu[i] - mu[i+test.md];
    field[i+test.md] = mu[i] + mu[i+test.md];
  }

  test.solve_eqn(field);

  for(int i = 0; i<= test.n_step; i++)
    std::cout<<test.ptnfn(i)<<std::endl;
  free(field);
  free(mu);

  return 0;
}

