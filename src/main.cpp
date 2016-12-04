#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "iterator.hpp"

int main(int argc, char * argv[])
{
  if(argc < 7)
  {
    std::cout<<"./main n_step bw m alpha beta kappa tau chiN domain."<<std::endl;
    return 0;
  }
  /*    int n_step = atoi(argv[1]); */
  /*    int bw = atoi(argv[2]); */
  /*    int m = atoi(argv[3]); */
     double alpha = atof(argv[1]);
     double beta = atof(argv[2]);
     double kappa = atof(argv[3]);
     double tau = atof(argv[4]);
     double chiN = atof(argv[5]);
     double domain = atof(argv[6]);

    char filename[] = "./data/transoutput_512_2_2_30_15_6_5";
    char transfilename[] = "./data/transoutput_512_2_2_40_20_6_5";
    char pdfname[] = "./data/pdf_512_2_2_40_20_6_5";
    char pdftransname[] = "./data/pdftrans_512_2_2_40_20_6_5";

    void (Iterator::*fp)();
    void (Iterator::*fp2)();
    fp = &Iterator::update_field;
    fp2 = &Iterator::delta_mu;
    Iterator test(400,8,512,alpha,beta,kappa,tau,chiN,domain);
    Iterator * obp=&test;

    test.read_data(filename);

    /* Picard pc; */
    SteepD sd;
    Anderson ad;

    //pc.solve(obp,fp,test.field,test.md*2);


    sd.solve(obp,fp2,test.mu,test.dmu,test.md*2,50);
    test.tensor();
    test.save_data(transfilename);

    sd.solve(obp,fp2,test.mu,test.dmu,test.md*2,50);
    test.tensor();
    test.save_data(transfilename);

    sd.solve(obp,fp2,test.mu,test.dmu,test.md*2,50);
    test.tensor();
    test.save_data(transfilename);

    sd.solve(obp,fp2,test.mu,test.dmu,test.md*2,50);
    test.tensor();
    test.save_data(transfilename);

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

