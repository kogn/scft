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
        KPSolver();
        ~KPSolver();
        void save_data(std::string);
        double ptnfn(int s = 0);
        int n_step;

        double Q;
        double nA;
        void density(const double * field);
        double * phi;
    private:
        void solve_eqn(const double *);
        static int count;
        static double * hist_forward;
        static double * hist_backward;

        double kappa, tau, alpha, beta;
        fftw_complex gamma[2];

        double volume;
        double domain[DIM];

        static fftw_complex * matrix;
        static fftw_complex * matrix1;

        void init_data();
        void pdf();
        void tensor();

        double t, dt;
        double * S[6];
        double * f;

        void onestep(const double *);
        void laplace(fftw_complex );
        void gradient(fftw_complex );
        void constant(fftw_complex ,const double *);
};

#endif //__KPSOLVER_H__
