#ifndef __KPSOLVER_H__
#define __KPSOLVER_H__

#ifdef __cplusplus 
extern "C" {
#endif //__cplusplus 
#include <fftw3.h>
#ifdef __cplusplus 
}
#endif //__cplusplus 

#include "kptransform.h"


class KPSolver : public S2_Space_trans, public S2_trans
{
    public:
        KPSolver(const Config &);
        ~KPSolver();
        void save_data(std::string);
        double ptnfn(int s = 0);
        int n_step;
        double f;

        double Q;
        void init_data_forward();
        void init_data_backward();
        void density();
        void pdf();
        void solve_eqn_forward(const double *);
        void solve_eqn_backward(const double *);
        bool head_tail;

        double * phi;
        double * S[6];
        double * dist;
    private:
        static int count;
        static double * hist_forward;
        static double * hist_backward;

        fftw_complex gamma[2];

        double volume;
        double domain[DIM];

        static fftw_complex * matrix;
        static fftw_complex * matrix1;

        void tensor();

        double t, dt;
        double alpha, beta;

        void onestep(const double *);
        void laplace(fftw_complex );
        void gradient(fftw_complex );
        void constant(fftw_complex ,const double *);
};

#endif //__KPSOLVER_H__
