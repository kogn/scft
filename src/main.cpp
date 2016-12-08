#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include "iterator.hpp"
#include "Config.h"

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


    void (Iterator<Solver,Solver>::*fp)();
    void (Iterator<Solver,Solver>::*fp2)();
    fp = &Iterator<Solver,Solver>::update_field;
    fp2 = &Iterator<Solver,Solver>::delta_mu;
    Iterator<Solver,Solver> test(configSettings);
    Iterator<Solver,Solver> * obp=&test;

    test.read_mu(input_filedir+input_filename);

    Picard pc;
    SteepD sd(output_filedir+output_filename);
    Anderson ad;

    //pc.solve(obp,fp,test.field,test.md*2);


    //sd.read_data2("./data/SteepD_20", test.mu, test.m[0]*2, test.m[1]);
    //sd.read_data("./data/SteepD_20", test.mu, test.md*2);
    sd.solve(obp,fp2,test.mu,test.dmu,test.md*2,max_steps);
    test.A.save_data(output_filedir+param_filename);


    /* ad.solve(obp,fp,test.field,test.md*2,20); */

  return 0;
}

