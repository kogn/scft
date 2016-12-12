#ifndef __ITERATOR_H__
#define __ITERATOR_H__

#include "solver.h"
#include <cmath>
#include <sstream>
#include <string>
#include <fstream>

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
class Iterator
{
    public:
        void update_field();
        void delta_mu();
        Iterator(const Config&);
        ~Iterator();
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
void Iterator<TA,TB>::read_field(std::string s)
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
void Iterator<TA,TB>::read_mu(std::string s)
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

class Anderson 
{
    public:
        Anderson();
        template<class T>
            void solve(T * ob,void (T::*func) (),double* ,int, int);
    private:
        double alpha;
        bool stop();
        double eps;
        int mk;
};

class Picard
{
    public:
        Picard();
        template<class T>
            void solve(T * ob,void (T::*func) (), double *, int, int);
    private:
        double alpha;
        double eps;
};
class SteepD
{
    public:
        void read_data(std::string filename, double * data, int length);
        void read_data2(std::string filename, double * data, int length, int times);
        void save_data(std::string filename, double * data, int length);
        std::string output_fileprefix;
        SteepD(std::string s);
        template<class T>
            void solve(T * ob,void (T::*func) (),double*, double *, int,int);
    private:
        double steplength;
        double eps;
};

static std::string num2str(int i)
{
    std::stringstream ss;
    ss << i;
    return ss.str();
}

    template<class T>
void SteepD::solve(T * ob,void (T::*func) (),double*x, double *dx, int n, int max_steps=200)
{
    int n_iters = 0;
    double err;
    timer();
    do{
        n_iters ++;
        std::cout<<"The "<<n_iters<<"th step:" <<std::endl;
        (ob->*func)();

        err = 0;
#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:err)
        for(int i = 0; i<n; i++)
        {
            x[i] -= steplength*dx[i];
            /* err = std::max(fabs(dx[i]),err); */
            err += dx[i]*dx[i];
        }
        err = sqrt(err/n);
        std::cout<<"Time = "<< timer()/60.<<" min, " <<"error = "<< err <<std::endl;
        save_data(output_fileprefix+"SteepD_"+num2str(n_iters), x, n);
    }while(err > eps && n_iters < max_steps);
    return;
}

    template<class T>
void Picard::solve(T * ob, void (T::*func)(), double * x, int n, int max_steps=200)
{
    int n_iters = 0;
    double * xcp = (double *)malloc(sizeof(double)*n);
    double err;
    do{
        memcpy(xcp, x, sizeof(double)*n);
        (ob->*func)();
        err = 0.;
#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:err)
        for(int i=0; i<n; i++)
        {
            //err = fabs(x[i]-xcp[i])>err?fabs(x[i]-xcp[i]):err;
            /* err = std::max(fabs(x[i]-xcp[i]),err); */
            err += (x[i]-xcp[i])*(x[i]-xcp[i]);
            x[i] = xcp[i]*(1.-alpha) + x[i]*alpha;
        }
        err = sqrt(err/n);
        n_iters ++;
        std::cout<<"The "<<n_iters<<"th step, "<<"error = "<< err <<std::endl;
    }while(n_iters<=max_steps&&err>eps);
    free(xcp);
    return;
}

    template<class T>
void Anderson::solve(T * ob,void (T::*func) (),double* y ,int n, int max_steps=200)
{
    int n_iters = 0;
    double err;
    double ** ls_A = (double**)malloc(sizeof(double*)*mk);
    double ** x = (double**)malloc(sizeof(double*)*mk);
    double ** g = (double**)malloc(sizeof(double*)*mk);

    ls_A[0] = (double *)malloc(sizeof(double)*n*mk);
    double * ls_b = (double *)malloc(sizeof(double)*n);
    g[0] = (double *)malloc(sizeof(double)*n*mk);
    x[0] = (double *)malloc(sizeof(double)*n*mk);
    for(int i = 1; i<mk; i++)
    {
        ls_A[i] = ls_A[0] + n*i;
        g[i] = g[0] + n*i;
        x[i] = x[0] + n*i;
    }

    do{
        int n_mod = n_iters%mk;
        err = 0.;
        memcpy(x[n_mod],y,sizeof(double)*n);
        (ob->*func)();
        memcpy(g[n_mod],y,sizeof(double)*n);
        for(int i=0; i<n; i++)
        {
            //err = std::max(err,fabs(y[i]-x[n_mod][i]));
            err += pow(fabs(y[i]-x[n_mod][i]),2);
        }
        err = sqrt(err/n);
        if(n_iters < mk-1){
            for(int i=0; i<n; i++)
                y[i] = y[i]*alpha + x[n_mod][i]*(1.-alpha);
        }
        else{
            for(int i = 0; i<n; i++)
            {
                ls_b[i] = g[mk-1][i] - x[mk-1][i];
                for(int j = 0; j<mk-1; j++)
                {
                    ls_A[j][i] = (g[j+1][i]-x[j+1][i]) - (g[j][i]-x[j][i]);
                }
            }
            LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', n, mk-1, 1, ls_A[0], n, ls_b, n);
            for(int i = 0; i<n; i++)
            {
                y[i] = g[mk-1][i];
                for(int j = 0; j<mk-1; j++)
                {
                    y[i]-=ls_b[j]*(g[j+1][i] - g[j][i]);
                }
            }
        }
        n_iters ++;
        std::cout<<"The "<<n_iters<<"th step, "<<"error = "<< err <<std::endl;
    }while(err > eps&&n_iters<max_steps);

    free(x[0]);
    free(g[0]);
    free(ls_A[0]);
    free(ls_b);
    free(x);
    free(g);
    free(ls_A);
    return;
}


template<class TA, class TB>
Iterator<TA,TB>::Iterator(const Config & configSettings):
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
                mu[i*m[1]+j] = -2*cos(2*M_PI*i/m[0]) + 2*cos(2*M_PI*j/m[1]);
                mu[i*m[1]+j+md] = -mu[i*m[1]+j];
                field[i*m[1]+j] = mu[i*m[1]+j] - mu[i*m[1]+j+md];
                field[i*m[1]+j+md] = mu[i*m[1]+j] + mu[i*m[1]+j+md];
            }
    }
}

template<class TA, class TB>
Iterator<TA,TB>::~Iterator()
{
  free(mu);
  free(dmu);
  free(field);
}


template<class TA, class TB>
void Iterator<TA,TB>::update_field()
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
void Iterator<TA, TB>::delta_mu()
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
double Iterator<TA,TB>::energy()
{
  H = 0;
  double h = 0;
#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:h)
  for(int i = 0; i<md; i++)
  {
    h += -mu[i] + mu[i+md]*mu[i+md]/chiN;
  }
  h /= md;
  H = h- (log(A.Q)*A.nA + log(B.Q)*B.nB);
  return H;
}


template<class TA, class TB>
void Iterator<TA,TB>::print_info()
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

#endif //__ITERATOR_H__
