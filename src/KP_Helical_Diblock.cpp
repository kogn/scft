#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include "Config.h"
#include "solver.h"
#include "kpsolver.h"
#include "diblock.hpp"
#include "diblock_melts.hpp"
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

    typedef Diblock<KPSolver,Solver> TA;

    void (Diblock_melts<TA>::*fp)();
    void (Diblock_melts<TA>::*fp2)();
    fp = &Diblock_melts<TA>::update_field;
    fp2 = &Diblock_melts<TA>::delta_mu;
    Diblock_melts<TA> test(configSettings);
    Diblock_melts<TA> * obp=&test;

    //test.read_mu(input_filedir+input_filename);
    //test.read_field(input_filedir+input_filename);

    Picard pc;
    SteepD sd(output_filedir+output_filename);
    Anderson ad(output_filedir+output_filename);
    //pc.solve(obp,fp,test.field,test.md*2);


    //sd.solve(obp,fp2,test.mu,test.dmu,test.md*2,max_steps);
    ad.solve(obp,fp,test.field,test.md*2,max_steps);
    test.A.save_data(output_filedir+param_filename);

  return 0;
}

