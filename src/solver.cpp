#include <cstdlib>
#include <cmath>
#include "solver.h"
#include <iostream>
#include <fstream>
#include <string>

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
#include "timer.h"

#ifdef __cplusplus
}
#endif //__cplusplus

#ifndef NUM_THREADS
#define NUM_THREADS 8
#endif //NUM_THREADS

int Solver::count = 0;
fftw_complex * Solver::matrix = NULL;
fftw_complex * Solver::matrix1 = NULL;
double * Solver::hist_forward = NULL;
//double * Solver::hist_backward = NULL;
Solver::Solver(const Config & configSettings):
    Space_trans(configSettings),SO3_trans(configSettings),Data(configSettings)
{
    n_step = configSettings.Read<int>("Steps_on_chain");
    alpha = configSettings.Read<double>("alpha");
    beta = configSettings.Read<double>("beta");
    kappa = configSettings.Read<double>("kappa");
    tau = configSettings.Read<double>("tau");
    phi = (double *) malloc(sizeof(double)*md);
    S[0] = (double *) malloc(sizeof(double)*md*6);
    for(int i = 0; i<6; i++)
    {
        S[i] = S[0] + md*i;
    }
    f = (double *)malloc(sizeof(double)*md*n3);
    dt = 1./n_step;
    gamma[0][0] = 0.324396404020171225;
    gamma[0][1] = 0.134586272490806680;
    gamma[1][0] = 0.351207191959657661;
    gamma[1][1] = -0.269172544981613415;
    /* gamma[0][0] = 1./3; */
    /* gamma[0][1] = 0.; */
    /* gamma[1][0] = 1./3; */
    /* gamma[1][1] = 0; */
    if(count == 0){
        matrix = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*totalCoeffs_so3(bw));
        matrix1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*totalCoeffs_so3(bw));
        getmatrix(bw,alpha,beta,kappa,tau,gamma, dt,matrix,matrix1);
        count ++;
        hist_forward = (double *) malloc(sizeof(double)*md*n3*(n_step+1));
        //hist_backward= (double *) malloc(sizeof(double)*md*n3*(n_step+1));
    }

    domain[0] = configSettings.Read<double>("domain0");
    if(DIM==2){
        domain[1] = configSettings.Read<double>("domain1");
    }
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
    if(count == 1){
        free(hist_forward);
        //free(hist_backward);
        count --;
        fftw_free(matrix);
        fftw_free(matrix1);
    }
    free(phi);
    free(S[0]);
    free(f);
}

void Solver::onestep(const double * field)
{
  fftw_complex dt_gamma0_2 = {dt*gamma[0][0]/2,dt*gamma[0][1]/2};
  fftw_complex dt_gamma0 = {dt*gamma[0][0],dt*gamma[0][1]};

  fftw_complex dt_gamma1_2 = {dt*gamma[1][0]/2,dt*gamma[1][1]/2};
  fftw_complex dt_gamma1 = {dt*gamma[1][0],dt*gamma[1][1]};

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
  /* for(int i = 0; i<md*n3; i++){ */
  /*     //realdata[i][0] = sqrt(realdata[i][0]*realdata[i][0]+realdata[i][1]*realdata[i][1]); */
  /*     //realdata[i][0] = fabs(realdata[i][0]); */
  /*     //realdata[i][1] = 0.; */
  /* } */

  t += dt;
}

void Solver::constant(fftw_complex dt,const double * field)
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
            int index1 = (i1+m[1]/2)%m[1]-m[1]/2;
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
                        int index = (i*m[1]+i1)*n3+n*n*j+n*k+l;
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

void Solver::solve_eqn(const double * field)
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

void Solver::pdf()
{
    double tmp;
    tmp = dt/(8.*M_PI*M_PI*Q);
    double * func = f;

#pragma omp parallel for num_threads(NUM_THREADS)
  for(int i = 0; i<md; i++)
    for(int j = 0; j<n; j++)
      for(int k = 0; k<n; k++)
        for(int kk = 0; kk<n; kk++)
        {
          int index = kk+n*k+j*n*n+i*n3;
          int index2 = (n-kk)%n+((n/2+k)%n)*n+(n-j-1)*n*n+i*n3;
          func[index] = 3./8*(hist_forward[md*n3*n_step+index2]+hist_forward[md*n3*n_step+index]);
          func[index] += 7./6*(hist_forward[md*n3*(n_step-1)+index2]*hist_forward[md*n3*1+index]
              +hist_forward[md*n3*1+index2]*hist_forward[md*n3*(n_step-1)+index]);
          func[index] += 23./24*(hist_forward[md*n3*(n_step-2)+index2]*hist_forward[md*n3*2+index]
              +hist_forward[md*n3*2+index2]*hist_forward[md*n3*(n_step-2)+index]);
          for(int ii = 3; ii<n_step-2; ii++)
          {
            func[index] += hist_forward[md*n3*ii+index]*hist_forward[md*n3*(n_step-ii)+index2];
          }
          func[index] = func[index]*tmp;
        }
  return;
}

void Solver::density(const double * field)
{
    solve_eqn(field);
    Q = ptnfn();
    pdf();
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<md; i++)
    {
        double tmp=0;
        double tmp1 = 0;
        for(int j = 0; j<n; j++)
        {
            double tmp2=0;
            double tmp3 =0;
            for(int k = 0; k<n*n; k++)
            {
                tmp2 += f[k+j*n*n+i*n3];
                tmp3 += f[k+j*n*n+i*n3];
            }
            tmp2 *= weights[j];
            tmp3 *= weights[j];
            tmp += tmp2;
      tmp1 += tmp3;
    }
    phi[i] = tmp * 4.*bw*M_PI*M_PI/n3;
    phi[i] = tmp1 * 4.*bw*M_PI*M_PI/n3;
  }
  return;
}

void Solver::tensor()
{
#pragma omp parallel for num_threads(NUM_THREADS)
  for(int i = 0; i<md; i++)
  {
    double tmp[6];
    double tmp1[6];
    for(int a = 0; a<6; a++)
    {
      tmp[a] = 0;
      tmp1[a] = 0;
    }
    for(int j = 0; j<n; j++)
    {
      double wt = weights[j];
      double s1 = sin(M_PI*(2*j+1)/4./bw);
      double c1 = cos(M_PI*(2*j+1)/4./bw);
      for(int k = 0; k<n; k++)
      {
        double c2 = cos(2*M_PI*k/n);
        double s2 = sin(2*M_PI*k/n);
        double m1 = c2*s1;
        double m2 = s2*s1;
        double m3 = c1;
        for(int l = 0; l<n; l++){
          int index =l+ k*n+j*n*n+i*n3;

          tmp[0] += f[index]*(m1*m1-1./3)*wt;

          tmp[1] += f[index]*(m2*m2-1./3)*wt;

          tmp[2] += f[index]*(m3*m3-1./3)*wt;

          tmp[3] += f[index]*m1*m2*wt;

          tmp[4] += f[index]*m1*m3*wt;

          tmp[5] += f[index]*m2*m3*wt;
        }
      }
    }
    for(int a = 0; a<6; a++)
    {
      S[a][i] = tmp[a] * 8.*bw*M_PI*M_PI/n3;
    }
  }
  return;
}

double Solver::ptnfn(int s)
{
  double q = 0;
  double * ptr = hist_forward+md*n3*s;
  double * ptr2 = hist_forward+md*n3*(n_step-s);
  //double * ptr3 = hist_backward+md*n3*s;
  //double * ptr4 = hist_backward+md*n3*(n_step-s);
#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:q)
  for(int i = 0; i<md; i++)
  {
    double tmp1=0;
    for(int j = 0; j<n; j++)
    {
      double tmp2=0;
      for(int k = 0; k<n; k++)
        for(int kk = 0; kk<n; kk++)
        {
          tmp2 += ptr[kk+n*k+j*n*n+i*n3]*ptr2[(n-kk)%n+((n/2+k)%n)*n+(n-j-1)*n*n+i*n3];
          //+ ptr3[kk+n*k+j*n*n+i*n3]*ptr4[(n-kk)%n+(n/2+k)%n*n+(n-j-1)*n*n+i*n3];
        }
      tmp2 *= weights[j];
      tmp1 += tmp2;
    }
    q += tmp1;
  }
  q = q*bw/(md*n3);
  return q;
}

void Solver::save_data(std::string filename)
{
  std::ofstream file(filename.c_str());
  tensor();
  for(int i = 0; i<md; i++)
  {
    file<<phi[i];
    for(int a = 0; a<6; a++)
      file<<" "<<S[a][i];
    file<<std::endl;
  }
  file.close();
  return;
}
