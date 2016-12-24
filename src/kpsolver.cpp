#include <cmath>
#include <iostream>

#include "kpsolver.h"

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus
#include <fftw3.h>
#include <lapacke.h>
#include "matrix.h"
#ifdef MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

#include <s2kit/makeweights.h>
#include <s2kit/cospmls.h>
#include <s2kit/FST_semi_memo.h>
#include <s2kit/csecond.h>

#include <soft/utils_so3.h>

#include <omp.h>
#ifdef __cplusplus
}
#endif //__cplusplus

int KPSolver::count = 0;
fftw_complex * KPSolver::matrix = NULL;
fftw_complex * KPSolver::matrix1 = NULL;
double * KPSolver::hist_forward = NULL;
double * KPSolver::hist_backward = NULL;
KPSolver::KPSolver(const Config & configSettings):
    S2_Space_trans(configSettings),S2_trans(configSettings),S2Data(configSettings)
{
    f = configSettings.Read<double>("f_semiflexible");
    n_step = configSettings.Read<int>("Steps_on_chain_semiflexible");
    alpha = configSettings.Read<double>("alpha");
    beta = configSettings.Read<double>("beta");
    head_tail = configSettings.Read<int>("head_tail");

    phi = (double *) malloc(sizeof(double)*md);
    S[0] = (double *) malloc(sizeof(double)*md*6);
    for(int i = 0; i<6; i++)
    {
        S[i] = S[0] + md*i;
    }
    dist = (double *)malloc(sizeof(double)*md*n2);
    dt = 1./n_step;
    gamma[0][0] = 0.324396404020171225;
    gamma[0][1] = 0.134586272490806680;
    gamma[1][0] = 0.351207191959657661;
    gamma[1][1] = -0.269172544981613415;
    if(count == 0){
        matrix = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*totalCoeffs_so3(bw));
        matrix1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*totalCoeffs_so3(bw));
        getmatrix(bw,alpha,beta,0.,0.,gamma, dt,matrix,matrix1);
        count ++;
        hist_forward = (double *) malloc(sizeof(double)*md*n2*(n_step+1));
        if(!head_tail){
            hist_backward= (double *) malloc(sizeof(double)*md*n2*(n_step+1));
        }
    }

    domain[0] = configSettings.Read<double>("domain0");
    if(DIM>=2){
        domain[1] = configSettings.Read<double>("domain1");
    }
    if(DIM>=3){
        domain[2] = configSettings.Read<double>("domain2");
    }
}

void KPSolver::init_data_forward()
{
    t = 0.;
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<n2*md; i++)
    {
        realdata0[i] = 1;
        realdata1[i] = 0;
        hist_forward[i] = 1.;
    }
    return;
}
void KPSolver::init_data_backward()
{
    t = 0.;
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<n2*md; i++)
    {
        realdata0[i] = 1;
        realdata1[i] = 0;
        hist_backward[i] = 1.;
    }
    return;
}


KPSolver::~KPSolver()
{
    if(count == 1){
        free(hist_forward);
        if(!head_tail){
            free(hist_backward);
        }
        fftw_free(matrix);
        fftw_free(matrix1);
    }
    count --;
    free(phi);
    free(S[0]);
    free(dist);
}

void KPSolver::onestep(const double * field)
{
  fftw_complex dt_gamma0_2 = {dt*gamma[0][0]/2,dt*gamma[0][1]/2};
  fftw_complex dt_gamma0 = {dt*gamma[0][0],dt*gamma[0][1]};

  fftw_complex dt_gamma1_2 = {dt*gamma[1][0]/2,dt*gamma[1][1]/2};
  fftw_complex dt_gamma1 = {dt*gamma[1][0],dt*gamma[1][1]};

  constant(dt_gamma0_2, field);
  for_space();
  gradient(dt_gamma0_2);
  for_s2();
  laplace(dt_gamma0);
  inv_s2();
  gradient(dt_gamma0_2);
  inv_space();
  constant(dt_gamma0_2, field);


  constant(dt_gamma1_2, field);
  for_space();
  gradient(dt_gamma1_2);
  for_s2();
  laplace(dt_gamma1);
  inv_s2();
  gradient(dt_gamma1_2);
  inv_space();
  constant(dt_gamma1_2, field);

  constant(dt_gamma0_2, field);
  for_space();
  gradient(dt_gamma0_2);
  for_s2();
  laplace(dt_gamma0);
  inv_s2();
  gradient(dt_gamma0_2);
  inv_space();
  constant(dt_gamma0_2, field);
  /* for(int i = 0; i<md*n2; i++){ */
  /*     //realdata[i][0] = sqrt(realdata[i][0]*realdata[i][0]+realdata[i][1]*realdata[i][1]); */
  /*     //realdata[i][0] = fabs(realdata[i][0]); */
  /*     //realdata[i][1] = 0.; */
  /* } */

  t += dt;
}


void KPSolver::constant(fftw_complex dt,const double * field)
{
#pragma omp parallel for num_threads(NUM_THREADS)
  for(int i = 0; i<md; i++)
  {
    double expw = exp(-field[i]*dt[0]);
    double cosw = cos(-field[i]*dt[1]);
    double sinw = sin(-field[i]*dt[1]);
    for(int j = 0; j<n2; j++)
    {
      double real = realdata0[i*n2+j];
      double imag = realdata1[i*n2+j];
      realdata0[i*n2+j] = expw*(cosw*real - sinw*imag);
      realdata1[i*n2+j] = expw*(cosw*imag + sinw*real);
    }
  }
  return;
}
#if DIM == 1
void KPSolver::gradient(fftw_complex dt)
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
                int index = i*n2+n*j+k;
                double real = realdata0[index];
                double imag = realdata1[index];
                realdata0[index] = (co*real - si*imag)*ex;
                realdata1[index] = (co*imag + si*real)*ex;
            }
        }
    }
    return;
}
#elif DIM == 2
void KPSolver::gradient(fftw_complex dt)
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
                    double tmp00 = -sin(M_PI*(2*j+1)/4./bw)*cos(2*M_PI*k/n)*2.*M_PI/domain[0]*dt[0];
                    double tmp01 = -sin(M_PI*(2*j+1)/4./bw)*cos(2*M_PI*k/n)*2.*M_PI/domain[0]*dt[1];
                    double tmp10 = -sin(M_PI*(2*j+1)/4./bw)*sin(2*M_PI*k/n)*2.*M_PI/domain[1]*dt[0];
                    double tmp11 = -sin(M_PI*(2*j+1)/4./bw)*sin(2*M_PI*k/n)*2.*M_PI/domain[1]*dt[1];
                    double co = cos(index0*tmp00+index1*tmp10);
                    double si = sin(index0*tmp00+index1*tmp10);
                    double ex = exp(-tmp01*index0-tmp11*index1);
                    int index = (i*m[1]+i1)*n2+n*j+k;
                    double real = realdata0[index];
                    double imag = realdata1[index];
                    realdata0[index] = (co*real - si*imag)*ex;
                    realdata1[index] = (co*imag + si*real)*ex;
                }
            }
        }
    }
    return;
}
#elif DIM == 3
void KPSolver::gradient(fftw_complex dt)
{
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i=0; i<m[0]; i++)
    {
        for(int i1 = 0; i1<m[1]; i1++){
            for(int i2 = 0; i2<m[2]; i2++){
                int index0 = (i+m[0]/2)%m[0]-m[0]/2;
                int index1 = (i1+m[1]/2)%m[1]-m[1]/2;
                int index2 = (i2+m[2]/2)%m[2]-m[2]/2;
                for(int j = 0; j<n; j++)
                {
                    for(int k = 0; k<n; k++)
                    {
                        double tmp00 = -sin(M_PI*(2*j+1)/4./bw)*cos(2*M_PI*k/n)*2.*M_PI/domain[0]*dt[0];
                        double tmp01 = -sin(M_PI*(2*j+1)/4./bw)*cos(2*M_PI*k/n)*2.*M_PI/domain[0]*dt[1];
                        double tmp10 = -sin(M_PI*(2*j+1)/4./bw)*sin(2*M_PI*k/n)*2.*M_PI/domain[1]*dt[0];
                        double tmp11 = -sin(M_PI*(2*j+1)/4./bw)*sin(2*M_PI*k/n)*2.*M_PI/domain[1]*dt[1];
                        double tmp20 = -cos(M_PI*(2*j+1)/4./bw)*2.*M_PI/domain[2]*dt[0];
                        double tmp21 = -cos(M_PI*(2*j+1)/4./bw)*2.*M_PI/domain[2]*dt[1];
                        double co = cos(index0*tmp00+index1*tmp10+index2*tmp20);
                        double si = sin(index0*tmp00+index1*tmp10+index2*tmp20);
                        double ex = exp(-tmp01*index0-tmp11*index1-tmp21*index2);
                        int index = (i*m[1]*m[2]+i1*m[2]+i2)*n2+n*j+k;
                        double real = realdata0[index];
                        double imag = realdata1[index];
                        realdata0[index] = (co*real - si*imag)*ex;
                        realdata1[index] = (co*imag + si*real)*ex;
                    }
                }
            }
        }
    }
    return;
}
#endif

void KPSolver::laplace(fftw_complex dt)
{
    fftw_complex * ptr;
    if(dt[1] > 0)
        ptr = matrix;
    else 
        ptr = matrix1;
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<md; i++)
        for(int l = 0; l<bw; l++)
            for(int j = -l; j<=l; j++)
            {
                int number = totalCoeffs_so3(l);
                int index = number+(j+l)*(2*l+1)+j+l;
                int index2 = n_coeff*i + seanindex(j,l,bw);
                spectdata0[index2]=ptr[index][0]*spectdata0[index2]
                    -ptr[index][1]*spectdata1[index2];
                spectdata1[index2]=ptr[index][0]*spectdata1[index2]
                    +ptr[index][1]*spectdata0[index2];
            }
    return;
}


void KPSolver::solve_eqn_forward(const double * field)
{
  for(int i = 1; i<n_step+1; i++)
  {
    onestep(field);
    cblas_dcopy(md*n2,realdata0,1,hist_forward+i*md*n2,1);
  }
  return;
}

void KPSolver::solve_eqn_backward(const double * field)
{
  for(int i = 1; i<n_step+1; i++)
  {
    onestep(field);
    cblas_dcopy(md*n2,realdata0,1,hist_backward+i*md*n2,1);
  }
  return;
}

void KPSolver::pdf()
{
    double tmp;
    tmp = dt*2./(4.*M_PI*Q);
    double * func = dist;
    if(head_tail){
        hist_backward = hist_forward;
    }

#pragma omp parallel for num_threads(NUM_THREADS)
  for(int i = 0; i<md; i++)
    for(int j = 0; j<n; j++)
      for(int k = 0; k<n; k++)
        {
          int index = k+j*n+i*n2;
          int index2 = ((n/2+k)%n)+(n-j-1)*n+i*n2;
          func[index] = 3./8*(hist_backward[md*n2*n_step+index2]+hist_forward[md*n2*n_step+index]);
          func[index] += 7./6*(hist_backward[md*n2*(n_step-1)+index2]*hist_forward[md*n2*1+index]
              +hist_backward[md*n2*1+index2]*hist_forward[md*n2*(n_step-1)+index]);
          func[index] += 23./24*(hist_backward[md*n2*(n_step-2)+index2]*hist_forward[md*n2*2+index]
              +hist_backward[md*n2*2+index2]*hist_forward[md*n2*(n_step-2)+index]);
          for(int ii = 3; ii<n_step-2; ii++)
          {
            func[index] += hist_forward[md*n2*ii+index]*hist_backward[md*n2*(n_step-ii)+index2];
          }
          func[index] = func[index]*tmp;
        }
  return;
}

void KPSolver::density()
{
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<md; i++)
    {
        double tmp=0;
        double tmp1 = 0;
        for(int j = 0; j<n; j++)
        {
            double tmp2=0;
            double tmp3 =0;
            for(int k = 0; k<n; k++)
            {
                tmp2 += dist[k+j*n+i*n2];
                tmp3 += dist[k+j*n+i*n2];
            }
            tmp2 *= weights[j];
            tmp3 *= weights[j];
            tmp += tmp2;
      tmp1 += tmp3;
    }
    phi[i] = tmp * 2.*bw*M_PI/n2;
    phi[i] = tmp1 * 2.*bw*M_PI/n2;
  }
  return;
}

void KPSolver::tensor()
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
            int index =k+j*n+i*n2;

            tmp[0] += dist[index]*(m1*m1-1./3)*wt;

            tmp[1] += dist[index]*(m2*m2-1./3)*wt;

            tmp[2] += dist[index]*(m3*m3-1./3)*wt;

            tmp[3] += dist[index]*m1*m2*wt;

            tmp[4] += dist[index]*m1*m3*wt;

            tmp[5] += dist[index]*m2*m3*wt;
        }
    }
    for(int a = 0; a<6; a++)
    {
        S[a][i] = tmp[a] *4.*bw*M_PI/n2;
    }
  }
  return;
}

double KPSolver::ptnfn(int s)
{
    double q = 0;
    double * ptr = hist_forward+md*n2*s;
    double * ptr2;
    if(head_tail){
        ptr2 = hist_forward+md*n2*(n_step-s);
    }else{
        ptr2 = hist_backward+md*n2*(n_step-s);
    }
    //double * ptr3 = hist_backward+md*n2*s;
  //double * ptr4 = hist_backward+md*n2*(n_step-s);
#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:q)
  for(int i = 0; i<md; i++)
  {
    double tmp1=0;
    for(int j = 0; j<n; j++)
    {
        double tmp2=0;
        for(int k = 0; k<n; k++)
            tmp2 += ptr[k+j*n+i*n2]*ptr2[((n/2+k)%n)+(n-j-1)*n+i*n2];
        //+ ptr3[kk+n*k+j*n*n+i*n2]*ptr4[(n-kk)%n+(n/2+k)%n*n+(n-j-1)*n*n+i*n2];
        tmp2 *= weights[j];
        tmp1 += tmp2;
    }
    q += tmp1;
  }
  Q = q*bw/(md*n2);
  return Q;
}

void KPSolver::save_data(std::string filename)
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

