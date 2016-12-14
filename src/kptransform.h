#ifndef __KPTRANSFORM_H__
#define __KPTRANSFORM_H__

#ifndef __cplusplus 
extern "C" {
#endif //__cplusplus
#include <fftw3.h>
#ifndef __cplusplus 
}
#endif //__cplusplus

#ifndef DIM
#define DIM 1
#endif
#include "Config.h"

class S2Data 
{
	public:
		int bw, n, n2, md, n_coeff;
        int m[DIM];
        double * realdata0;
        double * realdata1;
        double * spectdata0;
        double * spectdata1;
		S2Data(const Config &);
		virtual ~S2Data();
};

class S2_Space_trans:virtual public S2Data
{
	public:
		fftw_plan pf;
		fftw_plan pi;
		void for_space();
		void inv_space();

		S2_Space_trans(const Config &);
		~S2_Space_trans();
};

class S2_trans : virtual public S2Data
{
  public:
    void for_s2();
    void inv_s2();
    S2_trans(const Config &);
    ~S2_trans();
    double *weights ;
  private:
    int size, cutoff;
    double *seminaive_naive_tablespace, *trans_seminaive_naive_tablespace;
    double ** workspace;
    double **seminaive_naive_table, **trans_seminaive_naive_table;
    fftw_plan dctPlan, idctPlan;
    fftw_plan fftPlan, ifftPlan;
    fftw_iodim dims[1], howmany_dims[1];
    int rank, howmany_rank;
};
#endif //__KPTRANSFORM_H__
