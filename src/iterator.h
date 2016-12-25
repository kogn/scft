#ifndef __ITERATOR_H__
#define __ITERATOR_H__

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
class Anderson 
{
    public:
        void read_data(std::string filename, double * data, int length);
        void read_data2(std::string filename, double * data, int length, int times);
        void save_data(std::string filename, double * data, int length);
        std::string output_fileprefix;
        Anderson(std::string s);
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
    timer();

    do{
        std::cout<<"The "<<n_iters<<"th step:" <<std::endl;
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
        std::cout<<"Time = "<< timer()/60.<<" min, " <<"error = "<< err <<std::endl;
        save_data(output_fileprefix+"Anderson_"+num2str(n_iters), y, n);
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
#endif //__ITERATOR_H__
