#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include "Config.h"
#include "particle.h"
#include "homogeneous.hpp"
#include "solver.h"
#include "kpsolver.h"
#include "blends.hpp"
#include "iterator.h"

#ifndef DIM
#define DIM 1
#endif

int main(int argc, char * argv[])
{
    if(argc < 2)
    {
        std::cout<<"USE: ./main <path to config file>."<<std::endl;
        return 0;
    }
    Config configSettings(argv[1]);

    std::string output_filedir = configSettings.Read<std::string>("Output_dir");
    std::string input_filedir = configSettings.Read<std::string>("Input_dir");
    std::string input_filename = configSettings.Read<std::string>("Input_filename");
    std::string output_filename = configSettings.Read<std::string>("Output_filename");
    std::string param_filename= configSettings.Read<std::string>("Param_filename");
    int max_steps = configSettings.Read<int>("Max_steps");

    typedef Homogeneous<Solver> TA;
    typedef Particle TB;


    void (Iterator<TA,TB>::*fp)();
    void (Iterator<TA,TB>::*fp2)();
    fp = &Iterator<TA,TB>::update_field;
    fp2 = &Iterator<TA,TB>::delta_mu;
    Iterator<TA,TB> test(configSettings);
    Iterator<TA,TB> * obp=&test;

    //test.read_mu(input_filedir+input_filename);
    test.read_field(input_filedir+input_filename);

    Picard pc;
    SteepD sd(output_filedir+output_filename);
    Anderson ad(output_filedir+output_filename);
    //pc.solve(obp,fp,test.field,test.md*2);


    //sd.solve(obp,fp2,test.mu,test.dmu,test.md*2,max_steps);
    ad.solve(obp,fp,test.field,test.md*2,max_steps);
    test.A.save_data(output_filedir+param_filename);

  return 0;
}

