#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "kpsolver.h"
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


    double sum = 0;

    KPSolver test_solver(configSettings);

    for(int j = 0; j<test_solver.md; j++)
        for(int i = 0; i<test_solver.n_coeff; i++)
        {
            test_solver.spectdata0[test_solver.n_coeff*j+i] = j;
            test_solver.spectdata1[test_solver.n_coeff*j+i] = i;
        }
    test_solver.inv_s2();
    test_solver.inv_space();
    test_solver.for_space();
    test_solver.for_s2();

    test_solver.inv_s2();
    test_solver.inv_space();
    test_solver.for_space();
    test_solver.for_s2();

    test_solver.inv_s2();
    test_solver.inv_space();
    test_solver.for_space();
    test_solver.for_s2();

    for(int j = 0; j<test_solver.md; j++)
        for(int i = 0; i<test_solver.n_coeff; i++)
        {
            sum += fabs(test_solver.spectdata0[test_solver.n_coeff*j+i] - j);
            sum += fabs(test_solver.spectdata1[test_solver.n_coeff*j+i] - i);
        }
    std::cout<<sum<<std::endl;

    return 0;
}

