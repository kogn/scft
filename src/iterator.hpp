#ifndef __ITERATOR_H__
#define __ITERATOR_H__

#include "solver.h"
#include <cmath>

#ifdef __cplusplus 
extern "C" {
#endif //__cplusplus
#include <lapacke.h>
#include <memory.h>
#ifdef __cplusplus 
}
#endif //__cplusplus

#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif //NUM_THREADS

class Iterator : public Solver
{
	public:
		Iterator(int, int, int[], double, double, double, double,double,double[]);
		~Iterator();

		void update_field();
        void delta_mu();

		void print_info();

        void init();

		void quality();

        double energy();
		double H;

		double * mu;
        double * field;
        double * dmu;

		double chiN;

		void save_data(std::string);
		void save_pdf(std::string);
        void read_pdf(std::string);
		void read_data(std::string);

        void density();
        void pdf(double*);
        void tensor();

        double ptnfn(int);

        double Q_A, Q_B;
        double * S_A[6];
        double * S_B[6];
        double * phi_A, *phi_B;
        double * f_A, *f_B;

};

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
    SteepD();
    template<class T>
      void solve(T * ob,void (T::*func) (),double*, double *, int,int);
  private:
    double steplength;
    double eps;
};

template<class T>
void SteepD::solve(T * ob,void (T::*func) (),double*x, double *dx, int n, int max_steps=200)
{
  int n_iters = 0;
  double err;
  do{
    (ob->*func)();

    err = 0;
    for(int i = 0; i<n; i++)
    {
      x[i] -= steplength*dx[i];
      /* err = std::max(fabs(dx[i]),err); */
      err += dx[i]*dx[i];
    }
    err = sqrt(err/n);
    n_iters ++;
    std::cout<<"The "<<n_iters<<"th step, "<<"error = "<< err <<std::endl;
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

#endif //__ITERATOR_H__
