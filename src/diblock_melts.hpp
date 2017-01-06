#ifndef __DIBLOCK_MELTS__
#define __DIBLOCK_MELTS__

#include <cmath>
#include <sstream>
#include <string>
#include <fstream>
#include "Config.h"

#ifdef __cplusplus 
extern "C" {
#endif //__cplusplus
#include <lapacke.h>
#include <memory.h>
#include "timer.h"
#ifdef __cplusplus 
}
#endif //__cplusplus

#ifndef NUM_THREADS
#define NUM_THREADS 8
#endif //NUM_THREADS
#ifndef DIM 
#define DIM 1
#endif
template<class TA>
class Diblock_melts
{
    public:
        void update_field();
        void delta_mu();
        Diblock_melts(const Config&);
        ~Diblock_melts();
        double * mu;
        double * field;
        double * dmu;
        int md;
        int m[DIM];
        TA A;

        void read_field(std::string);
        void read_mu(std::string);
    private:
        void print_info();

        double energy();
        double H;
        double chiN;

};

    template<class TA>
void Diblock_melts<TA>::read_field(std::string s)
{
    std::ifstream file(s.c_str());
    for(int i = 0; i<md*2; i++)
    {
        file>>field[i];
    }
    for(int i = 0; i<md; i++){
        mu[i] = (field[i]+field[i+md])*0.5;
        mu[i+md] = (field[i+md]-field[i])*0.5;
    }
    file.close();
    return;
}
    template<class TA>
void Diblock_melts<TA>::read_mu(std::string s)
{
    std::ifstream file(s.c_str());
    for(int i = 0; i<md*2; i++)
    {
        file>>mu[i];
    }
    for(int i = 0; i<md; i++){
        field[i] = mu[i] - mu[i+md];
        field[i+md] = mu[i] + mu[i+md];
    }
    file.close();
    return;
}
template<class TA>
Diblock_melts<TA>::Diblock_melts(const Config & configSettings):
    A(configSettings)
{
    chiN = configSettings.Read<double>("chiN");
    md = A.md;

    mu= (double *)malloc(sizeof(double)*md*2);
    dmu= (double *)malloc(sizeof(double)*md*2);
    field = (double *)malloc(sizeof(double)*md*2);

    if(DIM == 1)
    {
        m[0] = A.m[0];
        for(int i = 0; i<md; i++)
        {
            mu[i] = -2*cos(2*M_PI*i/md);
            mu[i+md] = -mu[i];
            field[i] = mu[i] - mu[i+md];
            field[i+md] = mu[i] + mu[i+md];
        }
    }
    if(DIM == 2)
    {
        m[0] = A.m[0]; m[1] =  A.m[1];
        for(int i = 0; i<m[0]; i++)
            for(int j = 0; j<m[1]; j++)
            {
                /* field[i] = -2*cos(2*M_PI*i/md); */
                /* field[i+md] = -field[i]; */
                /* mu[i] = (field[i] + field[i+md])/2.; */
                /* mu[i+md] = (-field[i] + field[i+md])/2.; */
                /* mu[i] = -2*cos(2*M_PI*i/md); */
                /* mu[i+md] = -mu[i]; */
                /* field[i] = mu[i] - mu[i+md]; */
                /* field[i+md] = mu[i] + mu[i+md]; */
                mu[i*m[1]+j] = -.2*cos(2*M_PI*i/m[0]) + 2*cos(2*M_PI*j/m[1]);
                mu[i*m[1]+j+md] = -mu[i*m[1]+j];
                field[i*m[1]+j] = mu[i*m[1]+j] - mu[i*m[1]+j+md];
                field[i*m[1]+j+md] = mu[i*m[1]+j] + mu[i*m[1]+j+md];
            }
    }
    if(DIM == 3)
    {
        m[0] = A.m[0]; m[1] =  A.m[1]; m[2] =  A.m[2];
        for(int i = 0; i<m[0]; i++)
            for(int j = 0; j<m[1]; j++)
                for(int k = 0; k<m[2]; k++)
                {
                    mu[i*m[1]*m[2]+j*m[2]+k] = -.2*cos(2*M_PI*i/m[0]) + 2*cos(2*M_PI*j/m[1]) + .2*cos(2*M_PI*k/m[2]);
                    mu[i*m[1]*m[2]+j*m[2]+k+md] = -mu[i*m[1]*m[2]+j*m[2]+k];
                    field[i*m[1]*m[2]+j*m[2]+k] = mu[i*m[1]*m[2]+j*m[2]+k] - mu[i*m[1]*m[2]+j*m[2]+k+md];
                    field[i*m[1]*m[2]+j*m[2]+k+md] = mu[i*m[1]*m[2]+j*m[2]+k] + mu[i*m[1]*m[2]+j*m[2]+k+md];
                }
    }
}

    template<class TA>
Diblock_melts<TA>::~Diblock_melts()
{
    free(mu);
    free(dmu);
    free(field);
}


    template<class TA>
void Diblock_melts<TA>::update_field()
{
    A.density(field);
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0;i <md; i++)
    {
        double tmp;
        tmp = (field[i] + field[i+md] - chiN)*0.5;
        field[i] = chiN*A.B.phi[i]+tmp;
        field[i+md] = chiN*A.A.phi[i]+tmp;
        mu[i] = (field[i]+field[i+md])*0.5;
        mu[i+md] = (field[i+md]-field[i])*0.5;
    }
    energy();
    print_info();
    return;
}
    template<class TA>
void Diblock_melts<TA>::delta_mu()
{
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<md; i++)
    {
        field[i] = mu[i] - mu[i+md];
        field[i+md] = mu[i] + mu[i+md];
    }
    A.density(field);
    double tmp;
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<md; i++)
    {
        dmu[i] = -(A.A.phi[i]+A.B.phi[i]-1.);
        dmu[i+md]= A.B.phi[i] - A.A.phi[i] + 2*mu[i+md]/chiN;
    }
    energy();
    print_info();
    return;
}
    template<class TA>
double Diblock_melts<TA>::energy()
{
    H = 0;
    double h = 0;
#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:h)
    for(int i = 0; i<md; i++)
    {
        h += -mu[i] + mu[i+md]*mu[i+md]/chiN;
    }
    h /= md;
    H = h - log(A.Q);
    return H;
}


    template<class TA>
void Diblock_melts<TA>::print_info()
{
    double quality_A= 0.;
    double quality_B= 0.;
//#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:quality)
    for(int i = 0; i<A.md; i++)
    {
        quality_A += A.A.phi[i];
        quality_B += A.B.phi[i];
    }
    quality_A /= A.md;
    quality_B /= A.md;

    /* double maxQA = 0; */
    /* double minQA = 1000; */
    /* double sumQA = 0; */
    /* for(int i = 0; i<=A.A.n_step; i++){ */
    /*     //std::cout<<A.A.ptnfn(i)<<std::endl; */
    /*     A.A.ptnfn(i); */
    /*     sumQA += A.A.Q; */
    /*     maxQA = A.A.Q > maxQA ? A.A.Q:maxQA; */
    /*     minQA = A.A.Q < minQA? A.A.Q:minQA; */
    /* } */
    /* sumQA /= (A.A.n_step+1)*A.A.Q; */

    /* double maxQB = 0; */
    /* double minQB = 1e100; */
    /* double sumQB = 0; */
    /* for(int i = 0; i<=A.B.n_step; i++){ */
    /*     //std::cout<<A.B.ptnfn(i)<<std::endl; */
    /*     A.B.ptnfn(i); */
    /*     sumQB += A.B.Q; */
    /*     maxQB = A.B.Q > maxQB ? A.B.Q:maxQB; */
    /*     minQB = A.B.Q < minQB? A.B.Q:minQB; */
    /* } */
    /* sumQB /= (A.B.n_step+1)*A.B.Q; */
    std::cout<<"Quality = ("<<quality_A<<", "<<quality_B<<"); Q = "<<A.Q<<"; H = "
        <<H<<std::endl;

    /* std::cout<<"Quality = ("<<quality_A<<", "<<quality_B<<"); Q = ("<<A.A.ptnfn(0)<<", " */
    /*     <<A.A.ptnfn(A.A.n_step)<<", "<<A.B.ptnfn(0)<<", "<<A.B.ptnfn(A.B.n_step)<<"); H = " */
    return;
}

#endif //__DIBLOCK_MELTS__
