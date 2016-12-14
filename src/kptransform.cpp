#include <cmath>
#include <iostream>

#include "kptransform.h"


#ifdef __cplusplus
extern "C" {
#endif //__cplusplus
#include "matrix.h"

#include <lapacke.h>
#include <fftw3.h>
#include <s2kit/makeweights.h>
#include <s2kit/cospmls.h>
#include <s2kit/FST_semi_memo.h>
#include <s2kit/csecond.h>

#include <soft/utils_so3.h>
#include <omp.h>

#ifdef __cplusplus
}
#endif //__cplusplus

#ifndef NUM_THREADS
#define NUM_THREADS 1
#endif
#include "Config.h"


S2Data::S2Data(const Config & configSettings)
{
    bw = configSettings.Read<int>("Band_width");
    m[0] = configSettings.Read<int>("Grid_Size_x");
    if(DIM==2){
        m[1] = configSettings.Read<int>("Grid_Size_y");
    }
    n_coeff = bw*bw;
    n = bw*2;
    n2 = n*n;
    md = 1;
    for(int i = 0; i<DIM; i++)
        md *= m[i];
 
    realdata0 = (double *)malloc(sizeof(double)*n2*md);
    realdata1 = (double *)malloc(sizeof(double)*n2*md);
    spectdata0 = (double *)malloc(sizeof(double)*n_coeff*md);
    spectdata1 = (double *)malloc(sizeof(double)*n_coeff*md);
}

S2Data::~S2Data()
{
    free(realdata0);
    free(realdata1);
    free(spectdata0);
    free(spectdata1);
}

S2_Space_trans::S2_Space_trans(const Config & configSettings):S2Data(configSettings)
{
    fftw_iodim dims[DIM], howmany_dims[1];
    int rank = DIM;
    int howmany_rank = 1;
    if(DIM == 1){
        dims[0].n = m[0];
        dims[0].is = n2;
        dims[0].os = n2;
    }
    if(DIM == 2){
        dims[0].n = m[0];
        dims[0].is = n2*m[1];
        dims[0].os = n2*m[1];
        dims[1].n = m[1];
        dims[1].is = n2;
        dims[1].os = n2;
    }

    howmany_dims[0].n = n2;
    howmany_dims[0].is = 1;
    howmany_dims[0].os = 1;

    fftw_init_threads();
    fftw_plan_with_nthreads(NUM_THREADS);

    pf = fftw_plan_guru_split_dft( rank, dims,
            howmany_rank, howmany_dims,
            realdata0, realdata1,
            realdata0, realdata1,
            FFTW_ESTIMATE );

    pi = fftw_plan_guru_split_dft( rank, dims,
            howmany_rank, howmany_dims,
            realdata1, realdata0,
            realdata1, realdata0,
            FFTW_ESTIMATE );
    fftw_plan_with_nthreads(1);

    std::cout<<"space trans created!"<<std::endl;
}

S2_Space_trans::~S2_Space_trans()
{
    std::cout<<"space trans destoryed"<<std::endl;
    fftw_destroy_plan(pi);
    fftw_destroy_plan(pf);
}

void S2_Space_trans::for_space()
{
    fftw_execute(pf);
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<n2*md; i++)
    {
        realdata0[i]/=md;
        realdata1[i]/=md;
    }
    return;
}
void S2_Space_trans::inv_space()
{
    fftw_execute(pi);
    return;
}

S2_trans::S2_trans(const Config & configSettings):S2Data(configSettings)
{
    /*** ASSUMING WILL SEMINAIVE ALL ORDERS ***/
    cutoff = bw ;
    size = 2*bw;
    /* allocate lots of memory */
    seminaive_naive_tablespace =
        (double *) malloc(sizeof(double) *
                (Reduced_Naive_TableSize(bw,cutoff) +
                 Reduced_SpharmonicTableSize(bw,cutoff)));

    trans_seminaive_naive_tablespace =
        (double *) malloc(sizeof(double) *
                (Reduced_Naive_TableSize(bw,cutoff) +
                 Reduced_SpharmonicTableSize(bw,cutoff)));

    workspace = (double **)malloc(sizeof(double *)*NUM_THREADS);
    for(int i = 0; i<NUM_THREADS; i++)
    {
        workspace[i] = (double *) malloc(sizeof(double) * 
                ((8 * (bw*bw)) + 
                 (16 * bw)));
    }

    /* now precompute the Legendres */
    seminaive_naive_table = SemiNaive_Naive_Pml_Table(bw, cutoff,
            seminaive_naive_tablespace,
            workspace[0]);

    trans_seminaive_naive_table =
        Transpose_SemiNaive_Naive_Pml_Table(seminaive_naive_table,
                bw, cutoff,
                trans_seminaive_naive_tablespace,
                workspace[0]);

    /* make array for weights, and construct fftw plans */
    weights = (double *) malloc(sizeof(double) * 4 * bw);

    /* make DCT plans -> note that I will be using the GURU
       interface to execute these plans within the routines*/

    /* forward DCT */
    dctPlan = fftw_plan_r2r_1d( 2*bw, weights, realdata0,
            FFTW_REDFT10, FFTW_ESTIMATE ) ;

    /* inverse DCT */
    idctPlan = fftw_plan_r2r_1d( 2*bw, weights, realdata0,
            FFTW_REDFT01, FFTW_ESTIMATE );

    /*
       fftw "preamble" ;
       note that this plan places the output in a transposed array
       */
    rank = 1 ;
    dims[0].n = 2*bw ;
    dims[0].is = 1 ;
    dims[0].os = 2*bw ;
    howmany_rank = 1 ;
    howmany_dims[0].n = 2*bw ;
    howmany_dims[0].is = 2*bw ;
    howmany_dims[0].os = 1 ;
    /* forward fft */
    fftPlan = fftw_plan_guru_split_dft( rank, dims,
            howmany_rank, howmany_dims,
            realdata0, realdata1,
            workspace[0], workspace[0]+(4*bw*bw),
            FFTW_ESTIMATE );

    /*
       now plan for inverse fft - note that this plans assumes
       that I'm working with a transposed array, e.g. the inputs
       for a length 2*bw transform are placed every 2*bw apart,
       the output will be consecutive entries in the array
       */
    rank = 1 ;
    dims[0].n = 2*bw ;
    dims[0].is = 2*bw ;
    dims[0].os = 1 ;
    howmany_rank = 1 ;
    howmany_dims[0].n = 2*bw ;
    howmany_dims[0].is = 1 ;
    howmany_dims[0].os = 2*bw ;

    /* inverse fft */
    ifftPlan = fftw_plan_guru_split_dft( rank, dims,
            howmany_rank, howmany_dims,
            realdata0, realdata1,
            workspace[0], workspace[0]+(4*bw*bw),
            FFTW_ESTIMATE );


    /* now make the weights */
    makeweights( bw, weights );
}
void S2_trans::inv_s2()
{
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<md; i++)
    {
        int thread_num = omp_get_thread_num();
        /* do the inverse spherical transform */
        InvFST_semi_memo(spectdata0+i*n_coeff,spectdata1+i*n_coeff,
                realdata0+i*n2, realdata1+i*n2,
                bw,
                trans_seminaive_naive_table,
                workspace[thread_num],
                0,
                cutoff,
                &idctPlan,
                &ifftPlan );

    }
    return;
}
void S2_trans::for_s2()
{
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<md; i++)
    {
        int thread_num = omp_get_thread_num();
        FST_semi_memo(realdata0+i*n2, realdata1+i*n2,
                spectdata0+i*n_coeff, spectdata1+i*n_coeff,
                bw,
                seminaive_naive_table,
                workspace[thread_num],
                0,
                cutoff,
                &dctPlan,
                &fftPlan,
                weights );
    }
    return;
}

S2_trans::~S2_trans()
{

    fftw_destroy_plan( ifftPlan );
    fftw_destroy_plan( fftPlan );
    fftw_destroy_plan( idctPlan );
    fftw_destroy_plan( dctPlan );

    free( weights );
    free(trans_seminaive_naive_table);
    free(seminaive_naive_table);
    free(workspace);
    free(trans_seminaive_naive_tablespace);
    free(seminaive_naive_tablespace);
}


