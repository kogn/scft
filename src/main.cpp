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

    std::string output_filedir = configSettings.Read<std::string>("Output_dir");
    std::string input_filedir = configSettings.Read<std::string>("Input_dir");


    void (Iterator::*fp)();
    void (Iterator::*fp2)();
    fp = &Iterator::update_field;
    fp2 = &Iterator::delta_mu;
    Iterator test(n,bw,m,alpha,beta,kappa,tau,chiN,domain);
    Iterator * obp=&test;

    //test.read_data(input_filename);

    Picard pc;
    SteepD sd(output_filedir);
    Anderson ad;

    //pc.solve(obp,fp,test.field,test.md*2);


    sd.solve(obp,fp2,test.mu,test.dmu,test.md*2,50);
    test.tensor();
    test.save_data(output_filedir);

    sd.solve(obp,fp2,test.mu,test.dmu,test.md*2,50);
    test.tensor();
    test.save_data(output_filedir);

    sd.solve(obp,fp2,test.mu,test.dmu,test.md*2,50);
    test.tensor();
    test.save_data(output_filedir);


    /* ad.solve(obp,fp,test.field,test.md*2,20); */
    /* test.tensor(); */
    /* test.save_data(transfilename); */

    /* sd.solve(obp,fp2,test.mu,test.dmu,test.md*2,20); */
    /* test.tensor(); */
    /* test.save_data(transfilename); */
    //test.save_pdf(pdfname);

    /* test.read_pdf(pdftransname); */
    /* test.density(); */
    /* test.tensor(); */
    /* test.save_data(transfilename); */


    //test.save_pdf(pdfname);

  return 0;
}

