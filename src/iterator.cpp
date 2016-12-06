#include <cmath>
#include <iostream>
#include <fstream>

#include "iterator.hpp"
#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif //NUM_THREADS

#ifndef DIM 
#define DIM 1
#endif



Iterator::Iterator(int ns, int bw1, int m1[], double alpha0, double beta0, double kappa0, 
    double tau0, double chiN0,double domain0[]):
  Solver(ns,bw1,m1,alpha0,beta0,kappa0,tau0,domain0),Data(bw1,m1),chiN(chiN0)
{
  mu= (double *)malloc(sizeof(double)*md*2);
  dmu= (double *)malloc(sizeof(double)*md*2);
  field = (double *)malloc(sizeof(double)*md*2);
  phi_A = (double *) malloc(sizeof(double)*md);
  phi_B = (double *) malloc(sizeof(double)*md);
  S_A[0] = (double *) malloc(sizeof(double)*md*6);
  S_B[0] = (double *) malloc(sizeof(double)*md*6);
  for(int i = 0; i<6; i++)
  {
    S_A[i] = S_A[0] + md*i;
    S_B[i] = S_B[0] + md*i;
  }
  f_A = (double *)malloc(sizeof(double)*md*n3);
  f_B = (double *)malloc(sizeof(double)*md*n3);

  if(DIM == 1)
  {
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
    for(int i = 0; i<m[0]; i++)
      for(int j = 0; j<m[1]; j++)
      {
        mu[i*m[1]+j] = -2*cos(2*M_PI*i/m[0]) + 0.2*cos(2*M_PI*j/m[1]);
        mu[i*m[1]+j+md] = -mu[i*m[1]+j];
        field[i*m[1]+j] = mu[i*m[1]+j] - mu[i*m[1]+j+md];
        field[i*m[1]+j+md] = mu[i*m[1]+j] + mu[i*m[1]+j+md];
      }
  }
}

Iterator::~Iterator()
{
  free(mu);
  free(dmu);
  free(field);
  free(phi_A);
  free(phi_B);
  free(S_A[0]);
  free(S_B[0]);
  free(f_A);
  free(f_B);
}

double Iterator::ptnfn(int s=0)
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
          tmp2 += ptr[kk+n*k+j*n*n+i*n3]*ptr2[(n-kk)%n+(n/2+k)%n*n+(n-j-1)*n*n+i*n3];
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
void Iterator::pdf(double * func)
{
  double tmp;
  if(func == f_A)
    tmp = dt/(8.*M_PI*M_PI*Q_A);
  else 
    tmp = dt/(8.*M_PI*M_PI*Q_B);

#pragma omp parallel for num_threads(NUM_THREADS)
  for(int i = 0; i<md; i++)
    for(int j = 0; j<n; j++)
      for(int k = 0; k<n; k++)
        for(int kk = 0; kk<n; kk++)
        {
          int index = kk+n*k+j*n*n+i*n3;
          int index2 = (n-kk)%n+(n/2+k)%n*n+(n-j-1)*n*n+i*n3;
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

void Iterator::density()
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
      for(int k = 0; k<n*n; k++)
      {
        tmp2 += f_A[k+j*n*n+i*n3];
        tmp3 += f_B[k+j*n*n+i*n3];
      }
      tmp2 *= weights[j];
      tmp3 *= weights[j];
      tmp += tmp2;
      tmp1 += tmp3;
    }
    phi_A[i] = tmp * 4.*bw*M_PI*M_PI/n3;
    phi_B[i] = tmp1 * 4.*bw*M_PI*M_PI/n3;
  }
  return;
}

void Iterator::tensor()
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

          tmp[0] += f_A[index]*(m1*m1-1./3)*wt;
          tmp1[0]+= f_B[index]*(m1*m1-1./3)*wt;

          tmp[1] += f_A[index]*(m2*m2-1./3)*wt;
          tmp1[1]+= f_B[index]*(m2*m2-1./3)*wt;

          tmp[2] += f_A[index]*(m3*m3-1./3)*wt;
          tmp1[2]+= f_B[index]*(m3*m3-1./3)*wt;

          tmp[3] += f_A[index]*m1*m2*wt;
          tmp1[3]+= f_B[index]*m1*m2*wt;

          tmp[4] += f_A[index]*m1*m3*wt;
          tmp1[4]+= f_B[index]*m1*m3*wt;

          tmp[5] += f_A[index]*m2*m3*wt;
          tmp1[5]+= f_B[index]*m2*m3*wt;
        }
      }
    }
    for(int a = 0; a<6; a++)
    {
      S_A[a][i] = tmp[a] * 8.*bw*M_PI*M_PI/n3;
      S_B[a][i] = tmp1[a] *8.*bw*M_PI*M_PI/n3;
    }
  }
  return;
}

double Iterator::energy()
{
  H = 0;
  for(int i = 0; i<md; i++)
  {
    H += -mu[i] + mu[i+md]*mu[i+md]/chiN;
  }
  H /= md;
  H -= (log(Q_A) + log(Q_B))/2;
  return H;
}

void Iterator::update_field()
{
  solve_eqn(field);
  Q_A = ptnfn();
  pdf(f_A);

  solve_eqn(field+md);
  Q_B = ptnfn();
  pdf(f_B);

  density();
  for(int i = 0;i <md; i++)
  {
    double tmp;
    tmp = (field[i] + field[i+md] - chiN)*0.5;
    field[i] = chiN*phi_B[i]+tmp;
    field[i+md] = chiN*phi_A[i]+tmp;
    mu[i] = (field[i]+field[i+md])*0.5;
    mu[i+md] = (field[i+md]-field[i])*0.5;
  }
  energy();
  print_info();
  return;
}
void Iterator::delta_mu()
{
  for(int i = 0; i<md; i++)
  {
    field[i] = mu[i] - mu[i+md];
    field[i+md] = mu[i] + mu[i+md];
  }
  solve_eqn(field);
  Q_A = ptnfn();
  pdf(f_A);

  solve_eqn(field+md);
  Q_B = ptnfn();
  pdf(f_B);

  density();
  double tmp;
  for(int i = 0; i<md; i++)
  {
    dmu[i] = -(phi_A[i]+phi_B[i]-1.);
    dmu[i+md]= phi_B[i] - phi_A[i] +2*mu[i+md]/chiN;
  }
  energy();
  print_info();
  return;
}

void Iterator::save_pdf(std::string filename){
  std::ofstream file(filename.c_str());
  for(int i = 0; i<=md; i++)
    for(int j = -1; j<=n; j++)
      for(int k = 0; k<=n; k++)
        for(int kk = 0; kk<=n; kk++)
        {
          int index;
          if( j == -1)
            index = ((kk+n/2)%n)+n*((k+n/2)%n)+(i%md)*n3;
          else if(j == n)
            index = ((kk+n/2)%n)+n*((k+n/2)%n)+(n-1)*n*n+(i%md)*n3;
          else
            index = (kk%n)+n*(k%n)+(j%n)*n*n+(i%md)*n3;
          file<<f_B[index]<<std::endl;
        }
  file.close();
  return;
}

void Iterator::read_pdf(std::string filename){
  std::ifstream file(filename.c_str());
  for(int i = 0; i<md; i++)
    for(int j = 0; j<n; j++)
      for(int k = 0; k<n; k++)
        for(int kk = 0; kk<n; kk++)
        {
          int index = kk+n*k+j*n*n+i*n3;
          file>>f_B[index];
        }
  file.close();
  return;
}
void Iterator::save_data(std::string filename)
{
  std::ofstream file(filename.c_str());
  for(int i = 0; i<md; i++)
  {
    file<<phi_A[i]<<" "<<phi_B[i]<<" "<<field[i]<<" "<<field[i+md];
    for(int a = 0; a<6; a++)
      file<<" "<<S_A[a][i];
    for(int a = 0; a<6; a++)
      file<<" "<<S_B[a][i];
    file<<std::endl;
  }
  file.close();
  return;
}
void Iterator::read_data(std::string s)
{
  std::ifstream file(s.c_str());
  double tmp;
  for(int i = 0; i<md; i++)
  {
    file>>tmp>>tmp>>field[i]>>field[i+md];
    for(int j = 0; j<12; j++)
      file>>tmp;
    mu[i] = (field[i]+field[i+md])*0.5;
    mu[i+md] = (field[i+md]-field[i])*0.5;
  }
  file.close();
  return;
}

void Iterator::print_info()
{
    double quality= 0.;
    for(int i = 0; i<md; i++)
    {
        quality += phi_A[i] + phi_B[i];
    }
    quality /= md;
        std::cout<<std::endl<<"Quality = "<<quality<<"; Q_A = "<<Q_A<<
          "; Q_B = "<<Q_B<<"; H = "
        <<H<<std::endl;
    return;
}

Picard::Picard()
{
    alpha = 0.1;
    eps = 1e-4;
}

Anderson::Anderson()
{
    alpha = 0.1;
    eps = 1e-4;
    mk = 20;
}
SteepD::SteepD(std::string s):output_filedir(s)
{
  eps = 1e-4;
  steplength = 1;
}
