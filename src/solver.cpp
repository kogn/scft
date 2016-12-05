#include <cstdlib>
#include <cmath>
#include "solver.h"
#include <iostream>

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus
#include <fftw3.h>
#include <memory.h>
#include "matrix.h"
#include <cblas.h>

#include <soft/makeweights.h>
#include <soft/makeWigner.h>
#include <soft/utils_so3.h>
#include <soft/soft_fftw_pc.h>
#include <soft/csecond.h>
#include <omp.h>
#include <sys/time.h>

#ifdef __cplusplus
}
#endif //__cplusplus

#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif //NUM_THREADS

// Timer function: Measure the current time
double timer(void) {
    struct timeval Tvalue;
    struct timezone dummy;
    gettimeofday(&Tvalue, &dummy);
    double etime = (double)Tvalue.tv_sec + 1.0e-6*((double)Tvalue.tv_usec);
    return etime;
    //return omp_get_wtime();
}


Solver::Solver(int ns, int bw1, int m1[],double alpha0, double beta0, double kappa0, double tau0, double domain0[]):
  Space_trans(bw1,m1),SO3_trans(bw1,m1),Data(bw1,m1),n_step(ns),
  tau(tau0),alpha(alpha0),beta(beta0),kappa(kappa0)
{
  dt = 1./ns;
  gamma[0][0] = 0.324396404020171225;
  gamma[0][1] = 0.134586272490806680;
  gamma[1][0] = 0.351207191959657661;
  gamma[1][1] = -0.269172544981613415;
  /* gamma[0][0] = 1./3; */
  /* gamma[0][1] = 0.; */
  /* gamma[1][0] = 1./3; */
  /* gamma[1][1] = 0; */
  hist_forward = (double *) malloc(sizeof(double)*md*n3*(n_step+1));
  //hist_backward= (double *) malloc(sizeof(double)*md*n3*(n_step+1));
  matrix = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*totalCoeffs_so3(bw));
  matrix1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*totalCoeffs_so3(bw));
  getmatrix(bw,alpha,beta,kappa,tau,gamma, dt,matrix,matrix1);

  for(int i = 0; i<DIM; i++)
    domain[i] = domain0[i];
}

void Solver::init_data()
{
  t = 0.;
#pragma omp parallel for num_threads(NUM_THREADS)
  for(int i = 0; i<n3*md; i++)
  {
    realdata[i][0] = 1;
    realdata[i][1] = 0;
  }
  return;
}
Solver::~Solver()
{
  free(hist_forward);
  //free(hist_backward);
  fftw_free(matrix);
  fftw_free(matrix1);
}

void Solver::onestep(double * field)
{
  fftw_complex dt_gamma0_2 = {dt*gamma[0][0]/2,dt*gamma[0][1]/2};
  fftw_complex dt_gamma0 = {dt*gamma[0][0],dt*gamma[0][1]};

  fftw_complex dt_gamma1_2 = {dt*gamma[1][0]/2,dt*gamma[1][1]/2};
  fftw_complex dt_gamma1 = {dt*gamma[1][0],dt*gamma[1][1]};

  double etime, etime2;

  etime = timer();
  constant(dt_gamma0_2, field);
  for_space();
  gradient(dt_gamma0_2);
  for_so3();
  laplace(dt_gamma0);
  inv_so3();
  gradient(dt_gamma0_2);
  inv_space();
  constant(dt_gamma0_2, field);


  constant(dt_gamma1_2, field);
  for_space();
  gradient(dt_gamma1_2);
  for_so3();
  laplace(dt_gamma1);
  inv_so3();
  gradient(dt_gamma1_2);
  inv_space();
  constant(dt_gamma1_2, field);

  constant(dt_gamma0_2, field);
  for_space();
  gradient(dt_gamma0_2);
  for_so3();
  laplace(dt_gamma0);
  inv_so3();
  gradient(dt_gamma0_2);
  inv_space();
  constant(dt_gamma0_2, field);
  for(int i = 0; i<md; i++){
      realdata[i][0] = sqrt(realdata[i][0]*realdata[i][0]+realdata[i][1]*realdata[i][1]);
      realdata[i][1] = 0.;
  }

  t += dt;
}

void Solver::constant(fftw_complex dt,double * field)
{
#pragma omp parallel for num_threads(NUM_THREADS)
  for(int i = 0; i<md; i++)
  {
    double expw = exp(-field[i]*dt[0]);
    double cosw = cos(-field[i]*dt[1]);
    double sinw = sin(-field[i]*dt[1]);
    for(int j = 0; j<n3; j++)
    {
      double real = realdata[i*n3+j][0];
      double imag = realdata[i*n3+j][1];
      realdata[i*n3+j][0] = expw*(cosw*real - sinw*imag);
      realdata[i*n3+j][1] = expw*(cosw*imag + sinw*real);
    }
  }
  return;
}
#if DIM == 1
void Solver::gradient(fftw_complex dt)
{
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i=0; i<md; i++)
    {
        int index2 = (i+md/2)%md-md/2;
        for(int j = 0; j<n; j++)
        {
            double tmp = -cos(M_PI*(2*j+1)/4./bw)*2.*M_PI/domain[0]*dt[0];
            double tmp1 = -cos(M_PI*(2*j+1)/4./bw)*2.*M_PI/domain[0]*dt[1];
            double co = cos(index2*tmp);
            double si = sin(index2*tmp);
            double ex = exp(-tmp1*index2);
            for(int k = 0; k<n; k++)
            {
                //double tmp = -sin(M_PI*(2*j+1)/4./bw)*sin(2*M_PI*k/n)*2.*M_PI/domain[0]*dt[0];
                //double tmp1 = -sin(M_PI*(2*j+1)/4./bw)*sin(2*M_PI*k/n)*2.*M_PI/domain[0]*dt[1];
                //double tmp = -sin(M_PI*(2*j+1)/4./bw)*cos(2*M_PI*k/n)*2.*M_PI/domain[0]*dt[0];
                //double tmp1 = -sin(M_PI*(2*j+1)/4./bw)*cos(2*M_PI*k/n)*2.*M_PI/domain[0]*dt[1];
                for(int l = 0; l<n; l++)
                {
                    int index = i*n3+n*n*j+n*k+l;
                    double real = realdata[index][0];
                    double imag = realdata[index][1];
                    realdata[index][0] = (co*real - si*imag)*ex;
                    realdata[index][1] = (co*imag + si*real)*ex;
                }
            }
        }
    }
    return;
}
#elif DIM == 2
void Solver::gradient(fftw_complex dt)
{
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i=0; i<m[0]; i++)
    {
        for(int i1 = 0; i1<m[1]; i1++){
            int index0 = (i+m[0]/2)%m[0]-m[0]/2;
            int index1 = (i+m[1]/2)%m[1]-m[1]/2;
            for(int j = 0; j<n; j++)
            {
                for(int k = 0; k<n; k++)
                {
                    double tmp00 = -sin(M_PI*(2*j+1)/4./bw)*sin(2*M_PI*k/n)*2.*M_PI/domain[0]*dt[0];
                    double tmp01 = -sin(M_PI*(2*j+1)/4./bw)*sin(2*M_PI*k/n)*2.*M_PI/domain[0]*dt[1];
                    double tmp10 = -sin(M_PI*(2*j+1)/4./bw)*cos(2*M_PI*k/n)*2.*M_PI/domain[1]*dt[0];
                    double tmp11 = -sin(M_PI*(2*j+1)/4./bw)*cos(2*M_PI*k/n)*2.*M_PI/domain[1]*dt[1];
                    double co = cos(index0*tmp00+index1*tmp10);
                    double si = sin(index0*tmp00+index1*tmp10);
                    double ex = exp(-tmp01*index0-tmp11*index1);
                    for(int l = 0; l<n; l++)
                    {
                        int index = i*n3+n*n*j+n*k+l;
                        double real = realdata[index][0];
                        double imag = realdata[index][1];
                        realdata[index][0] = (co*real - si*imag)*ex;
                        realdata[index][1] = (co*imag + si*real)*ex;
                    }
                }
            }
        }
    }
    return;
}
#endif

void Solver::laplace(fftw_complex dt)
{
    fftw_complex * ptr;
    if(dt[1] > 0)
        ptr = matrix;
    else 
        ptr = matrix1;
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<md; i++)
    {
        for(int l = 0; l<bw; l++)
        {
            int number = totalCoeffs_so3(l);
            for(int j = -l; j<=l; j++)
                for(int k = -l; k<=l; k++)
                {
                    double real = 0.;
                    double imag = 0.;
                    for(int kk = -l;kk<=l;kk++)
                    {
                        int index = so3CoefLoc(j,kk,l,bw) + i*n_coeff;
                        int index2 = number +(k+l)*(2*l+1)+kk+l;
                        real+=ptr[index2][0]*spectdata[index][0]
                            -ptr[index2][1]*spectdata[index][1];
                        imag+=ptr[index2][0]*spectdata[index][1]
                            +ptr[index2][1]*spectdata[index][0];
                    }
                    spectdata[n_coeff*i+so3CoefLoc(j,k,l,bw)][0] = real;
                    spectdata[n_coeff*i+so3CoefLoc(j,k,l,bw)][1] = imag;
                }
        }
    }
    return;
}

void Solver::solve_eqn(double * field)
{
  init_data();

#pragma omp parallel for num_threads(NUM_THREADS)
  for(int j = 0; j<md*n3; j++)
    hist_forward[j] = 1.;

  for(int i = 1; i<n_step+1; i++)
  {
    onestep(field);
    cblas_dcopy(md*n3,realdata[0],2,hist_forward+i*md*n3,1);
/* #pragma omp parallel for num_threads(NUM_THREADS) */
/*     for(int j = 0; j<md*n3; j++) */
/*     { */
/*       hist_forward[j+i*md*n3] = realdata[j][0]; */
/*     } */
  }
  return;
}

