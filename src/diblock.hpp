#ifndef __DIBLOCK_HPP__
#define __DIBLOCK_HPP__
#include <cmath>

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus
#include <fftw3.h>

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
#include "kpsolver.h"
#include "solver.h"

template<typename TA, typename TB>
class Diblock
{
    public:
        TA A;
        TB B;
        Diblock(const Config &);
        //~Diblock();
        void save_data(std::string);
        double ptnfn(int s = 0);
        int n_step;

        double Q;
        double prop;
        int * m;
        int md;

        void density(const double * field);
        void solve_eqn(const double *);

    private:
        void pdf();
};

void connect_bond_forward(KPSolver & A, Solver & B){
    double * ptrA = A.hist_forward+A.md*A.n2*A.n_step;
    double * ptrB = B.hist_forward;
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<A.md*A.n2; i++)
    {
        for(int l = 0; l<A.n; l++){
            ptrA[i*A.n+l] = ptrA[i]/2./M_PI;
        }
    }
    return;
}

void connect_bond_backward(Solver & A, KPSolver & B){
    double * ptrB = B.hist_backward+B.md*B.n2*B.n_step;
    double * ptrA = A.hist_backward;
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<B.md*B.n2; i++)
    {
        for(int l = 0; l<B.n; l++){
            ptrA[i*A.n+l] = ptrA[i]/2./M_PI;
        }
    }
    return;
}

void connect_bond_forward(Solver & A, KPSolver & B){
    double * ptrA = A.hist_forward+A.md*A.n3*A.n_step;
    double * ptrB = B.hist_forward;
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<A.md*B.n2; i++)
    {
        ptrB[i] = 0.;
        for(int l = 0; l<A.n; l++){
            ptrB[i] += ptrA[i*A.n+l];
        }
        ptrB[i] *= 2.*M_PI/A.n;
    }
    return;
}


void connect_bond_backward(KPSolver& A, Solver & B){
    double * ptrB = B.hist_forward+B.md*B.n3*B.n_step;
    double * ptrA = A.hist_forward;
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<A.md*A.n2; i++)
    {
        ptrA[i] = 0.;
        for(int l = 0; l<A.n; l++){
            ptrA[i] += ptrB[i*A.n+l];
        }
        ptrA[i] *= 2.*M_PI/B.n;
    }
    return;
}

    template<typename TA, typename TB>
Diblock<TA,TB>::Diblock(const Config &configSettings):A(configSettings),B(configSettings)
{
    n_step = configSettings.Read<int>("Steps_on_chain");
    m = A.m;
    md = A.md;
}

    template<typename TA, typename TB>
void Diblock<TA,TB>::solve_eqn(const double * field)
{
    A.init_data_forward();
    A.solve_eqn_forward(field);
    connect_bond_forward(A,B);
    B.solve_eqn_forward(field+md);


    B.init_data_backward();
    B.solve_eqn_backward(field+md);
    connect_bond_backward(A,B);
    A.solve_eqn_backward(field);
    return;
}

    template<typename TA, typename TB>
void Diblock<TA,TB>::pdf()
{
    A.pdf();
    B.pdf();
    return;
}
    template<typename TA, typename TB>
void Diblock<TA,TB>::density(const double * field)
{
    solve_eqn(field);
    ptnfn();
    pdf();
    A.density();
    B.density();
    return;
}

    template<typename TA, typename TB>
double Diblock<TA,TB>::ptnfn(int s)
{
    if(s <= A.n_step){
        Q = A.ptnfn(s);
        B.Q = Q;
    }else{
        Q = B.ptnfn(s-A.n_step);
        A.Q = Q;
    }
    return Q;
}

    template<typename TA, typename TB>
void Diblock<TA,TB>::save_data(std::string filename)
{
    A.save_data(filename);
    B.save_data(filename);
    return;
}
#endif //__DIBLOCK_HPP__

