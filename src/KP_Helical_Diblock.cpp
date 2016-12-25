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

    std::string output_prefix= configSettings.Read<std::string>("Output_prefix");
    bool read_mu_from_file = configSettings.Read<int>("Read_mu_from_file");
    bool read_field_from_file = configSettings.Read<int>("Read_field_from_file");
    std::string input_filename = configSettings.Read<std::string>("Input_filename");
    std::string output_filename = configSettings.Read<std::string>("Output_filename");
    std::string param_filename= configSettings.Read<std::string>("Param_filename");
    int max_steps_SD = configSettings.Read<int>("Max_steps_SD");
    int max_steps_Ad = configSettings.Read<int>("Max_steps_Ad");

    typedef Diblock<Solver,KPSolver> TA;

    void (Diblock_melts<TA>::*fp)();
    void (Diblock_melts<TA>::*fp2)();
    fp = &Diblock_melts<TA>::update_field;
    fp2 = &Diblock_melts<TA>::delta_mu;
    Diblock_melts<TA> test(configSettings);
    Diblock_melts<TA> * obp=&test;

    if(read_mu_from_file){
        test.read_mu(input_filename);
    }
    if(read_field_from_file){
        test.read_field(input_filename);
    }

    //Picard pc;
    SteepD sd(output_prefix+output_filename);
    Anderson ad(output_prefix+output_filename);
    //pc.solve(obp,fp,test.field,test.md*2);


    if(max_steps_SD>0){
        sd.solve(obp,fp2,test.mu,test.dmu,test.md*2,max_steps_SD);
    }
    if(max_steps_Ad>0){
        ad.solve(obp,fp,test.field,test.md*2,max_steps_Ad);
    }
    test.A.save_data(output_prefix+param_filename);

  return 0;
}

