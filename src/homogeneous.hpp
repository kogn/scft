#ifndef __HOMOGENEOUS_HPP__
#define __HOMOGENEOUS_HPP__

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus
#include <fftw3.h>
#include <cblas.h>

#ifdef __cplusplus
}
#endif //__cplusplus
#include<string>

#ifndef NUM_THREADS
#define NUM_THREADS 8
#endif //NUM_THREADS
#ifndef DIM 
#define DIM  1
#endif

#include "transform.h"

template<typename TA>
class Homogeneous
{
    public:
        TA A;
        Homogeneous(const Config &);
        ~Homogeneous();
        void save_data(std::string);
        double ptnfn(int s = 0);
        int n_step;

        double Q;
        double prop;
        double * phi;
        int * m;
        int md;

        void density(const double * field);
        void solve_eqn(const double *);

    private:
        void pdf();
};
    template<typename TA>
Homogeneous<TA>::Homogeneous(const Config &configSettings):A(configSettings)
{
    n_step = configSettings.Read<int>("Steps_on_chain");
    prop = configSettings.Read<double>("nA");
    m = A.m;
    md = A.md;
    phi = new double[md];
}

template<typename TA>
Homogeneous<TA>::~Homogeneous(){
  delete [] phi;
}

    template<typename TA>
void Homogeneous<TA>::solve_eqn(const double * field)
{
    A.init_data_forward();
    A.solve_eqn_forward(field);
    return;
}

    template<typename TA>
void Homogeneous<TA>::pdf()
{
    A.pdf();
    return;
}
    template<typename TA>
void Homogeneous<TA>::density(const double * field)
{
    solve_eqn(field);
    ptnfn();
    pdf();
    A.density();
    //cblas_daxpy(md,prop,A.phi,1,phi,1);
    for(int i = 0; i<md; i++)
        phi[i] = prop * A.phi[i];
    return;
}

    template<typename TA>
double Homogeneous<TA>::ptnfn(int s)
{
    Q = A.ptnfn(s);
    return Q;
}

    template<typename TA>
void Homogeneous<TA>::save_data(std::string filename)
{
    A.save_data(filename);
    return;
}
#endif //__HOMOGENEOUS_HPP__
