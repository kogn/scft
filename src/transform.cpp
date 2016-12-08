#include <cmath>

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

#include <memory.h>
#include <omp.h>
#include <malloc.h>
#include <fftw3.h>
#include <lapacke.h>
#include <cblas.h>

#include <soft/makeweights.h>
#include <soft/makeWigner.h>
#include <soft/utils_so3.h>
#include <soft/soft_fftw_pc.h>
#include <soft/csecond.h>

#ifdef __cplusplus
}
#endif //__cplusplus

#ifndef DIM
#define DIM 1
#endif

#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif //NUM_THREADS

#include "transform.h"
#include "Config.h"

Data::Data(const Config & configSettings)
{
    bw = configSettings.Read<int>("Band_width");
    m[0] = configSettings.Read<int>("Grid_Size_x");
    if(DIM==2){
        m[1] = configSettings.Read<int>("Grid_Size_y");
    }

    n_coeff = totalCoeffs_so3(bw);
    n = bw*2;
    n3 = n*n*n;
    md = 1;
    for(int i = 0; i<DIM; i++)
        md *= m[i];
    realdata = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*n3*md);
    spectdata = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*n_coeff*md);
}

Data::~Data()
{
  fftw_free(realdata);
  fftw_free(spectdata);
}

Space_trans::Space_trans(const Config & configSettings):Data(configSettings)
{
  fftw_iodim dims[DIM], howmany_dims[1];
  int rank = DIM;
  int howmany_rank = 1;
  if(DIM == 1){
    dims[0].n = m[0];
    dims[0].is = n3;
    dims[0].os = n3;
  }
  if(DIM == 2){
    dims[0].n = m[0];
    dims[0].is = n3;
    dims[0].os = n3;
    dims[1].n = m[1];
    dims[1].is = n3*m[0];
    dims[1].os = n3*m[0];
  }
  howmany_dims[0].n = n3;
  howmany_dims[0].is = 1;
  howmany_dims[0].os = 1;

  fftw_init_threads();
  fftw_plan_with_nthreads(NUM_THREADS);

  pf = fftw_plan_guru_dft( rank, dims,
      howmany_rank, howmany_dims,
      realdata, realdata,
      FFTW_FORWARD, FFTW_ESTIMATE );

  pi = fftw_plan_guru_dft( rank, dims,
      howmany_rank, howmany_dims,
      realdata, realdata,
      FFTW_BACKWARD, FFTW_ESTIMATE );
  fftw_plan_with_nthreads(1);
}

Space_trans::~Space_trans()
{
  fftw_destroy_plan(pi);
  fftw_destroy_plan(pf);
}

void Space_trans::for_space()
{
  fftw_execute(pf);
  cblas_dscal(n3*md*2,1./md,realdata[0],1);

  /* for(int i = 0; i<n3*md; i++) */
  /* { */
  /*   realdata[i][0]/=md; */
  /*   realdata[i][1]/=md; */
  /* } */
  return;
}
void Space_trans::inv_space()
{
  fftw_execute(pi);
  return;
}

SO3_trans::SO3_trans(const Config & configSettings):Data(configSettings)
{
  workspace_cx =(fftw_complex**)malloc(sizeof(fftw_complex*)*NUM_THREADS);
  workspace_cx2=(fftw_complex**)malloc(sizeof(fftw_complex*)*NUM_THREADS);
  workspace_re = (double **)malloc(sizeof(double*)*NUM_THREADS);

  for(int i = 0; i<NUM_THREADS; i++)
  {
    workspace_cx[i] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n3);
    workspace_cx2[i] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n3);
    workspace_re[i] = (double *)malloc(sizeof(double)*(24*bw + 2*bw*bw));
  }

  wignerSpace = ((bw*bw)*(2+3*bw+bw*bw))/3 ;
  wigners = (double *)malloc(sizeof(double) * wignerSpace);
  wignersTrans = (double *)malloc(sizeof(double) * wignerSpace);

  weights = (double *) malloc(sizeof(double) * (2*bw));

  genWigAll( bw, wigners, workspace_re[0] );
  genWigAllTrans( bw, wignersTrans, workspace_re[0]);
  makeweights2( bw, weights );

  {
    int na[2], inembed[2], onembed[2];
    int rank, howmany, istride, idist, ostride, odist;
    howmany = n*n;
    idist = n;
    odist = n;
    rank = 2 ;
    inembed[0] = n;
    inembed[1] = n*n;
    onembed[0] = n;
    onembed[1] = n*n;
    istride = 1;
    ostride = 1;
    na[0] = 1;
    na[1] = n;

    p1 = fftw_plan_many_dft( rank, na, howmany,
        workspace_cx2[0], inembed,
        istride, idist,
        workspace_cx[0], onembed,
        ostride, odist,
        FFTW_BACKWARD, FFTW_MEASURE );
        //FFTW_BACKWARD, FFTW_PATIENT);

    p2 = fftw_plan_many_dft( rank, na, howmany,
        workspace_cx[0], inembed,
        istride, idist,
        workspace_cx2[0], onembed,
        ostride, odist,
        FFTW_FORWARD, FFTW_MEASURE );
        //FFTW_FORWARD, FFTW_PATIENT);
  }
}
void SO3_trans::inv_so3()
{
#pragma omp parallel for num_threads(NUM_THREADS)
  for(int i = 0; i<md; i++)
  {
    int thread_num = omp_get_thread_num();
    Inverse_SO3_Naive_fftw_pc(bw,
        spectdata+i*n_coeff,
        realdata+i*n3,
        workspace_cx[thread_num],
        workspace_cx2[thread_num],
        workspace_re[thread_num],
        &p2,
        wignersTrans,
        0 );
  }
  return;
}
void SO3_trans::for_so3()
{
#pragma omp parallel for num_threads(NUM_THREADS)
  for(int i = 0; i<md; i++)
  {
    int thread_num = omp_get_thread_num();
    Forward_SO3_Naive_fftw_pc( bw,
        realdata+i*n3,
        spectdata+i*n_coeff,
        workspace_cx[thread_num],
        workspace_cx2[thread_num],
        workspace_re[thread_num],
        weights,
        &p1,
        wigners,
        0 );
  }
  return;
}

SO3_trans::~SO3_trans()
{
  fftw_destroy_plan(p2);
  fftw_destroy_plan(p1);

  free(weights);
  free(wignersTrans) ;
  free(wigners) ;

  for(int i=0; i<NUM_THREADS; i++)
  {
    fftw_free(workspace_cx[i]);
    fftw_free(workspace_cx2[i]);
    free(workspace_re[i]);
  }
  free(workspace_re);
  free(workspace_cx2);
  free(workspace_cx);
}

