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
template<class TA, class TB>
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
        TB B;

        void read_field(std::string);
        void read_mu(std::string);
    private:
        void print_info();

        void quality();

        double energy();
        double H;
        double chiN;

};

    template<class TA, class TB>
void Diblock_melts<TA,TB>::read_field(std::string s)
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
    template<class TA, class TB>
void Diblock_melts<TA,TB>::read_mu(std::string s)
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
template<class TA, class TB>
Diblock_melts<TA,TB>::Diblock_melts(const Config & configSettings):
    A(configSettings),B(configSettings)
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

    template<class TA, class TB>
Diblock_melts<TA,TB>::~Diblock_melts()
{
    free(mu);
    free(dmu);
    free(field);
}


    template<class TA, class TB>
void Diblock_melts<TA,TB>::update_field()
{
    A.density(field);
    B.density(field+md);
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0;i <md; i++)
    {
        double tmp;
        tmp = (field[i] + field[i+md] - chiN)*0.5;
        field[i] = chiN*B.phi[i]+tmp;
        field[i+md] = chiN*A.phi[i]+tmp;
        mu[i] = (field[i]+field[i+md])*0.5;
        mu[i+md] = (field[i+md]-field[i])*0.5;
    }
    energy();
    print_info();
    return;
}
    template<class TA, class TB>
void Diblock_melts<TA, TB>::delta_mu()
{
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<md; i++)
    {
        field[i] = mu[i] - mu[i+md];
        field[i+md] = mu[i] + mu[i+md];
    }
    A.density(field);
    B.density(field+md);
    double tmp;
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<md; i++)
    {
        dmu[i] = -(A.phi[i]+B.phi[i]-1.);
        dmu[i+md]= B.phi[i] - A.phi[i] + 2*mu[i+md]/chiN;
    }
    energy();
    print_info();
    return;
}
    template<class TA, class TB>
double Diblock_melts<TA,TB>::energy()
{
    H = 0;
    double h = 0;
#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:h)
    for(int i = 0; i<md; i++)
    {
        h += -mu[i] + mu[i+md]*mu[i+md]/chiN;
    }
    h /= md;
    H = h- (log(A.Q)*A.prop + log(B.Q)*B.prop);
    return H;
}


    template<class TA, class TB>
void Diblock_melts<TA,TB>::print_info()
{
    double quality= 0.;
#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:quality)
    for(int i = 0; i<md; i++)
    {
        quality += A.phi[i] + B.phi[i];
    }
    quality /= md;
    std::cout<<"Quality = "<<quality<<"; Q_A = "<<A.Q<<
        "; Q_B = "<<B.Q<<"; H = "
        <<H<<std::endl;
    return;
}

#endif //__DIBLOCK_MELTS__
