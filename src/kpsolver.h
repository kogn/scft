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

/* template<typename TA, typename TB> */
/* void connect_bond_forward(TA&,TB&); */
/* template<typename TA, typename TB> */
/* void connect_bond_backward(TA&,TB&); */

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
        void pdfs();
        void solve_eqn_forward(const double *);
        void solve_eqn_backward(const double *);
        /* template<typename TA, typename TB> */
        /*     friend void connect_bond_forward<>(TA&,TB&); */
        /* template<typename TA, typename TB> */
        /*     friend void connect_bond_backward<>(TA&,TB&); */
        bool head_tail;

        double * phi;
        double * S[6];
        double * dist;
        double * dists;
        static double * hist_forward;
        static double * hist_backward;
    private:
        static int count;

        double lambda;
        fftw_complex gamma[2];

        double volume;
        double domain[DIM];


        void tensor();

        double t, dt;

        void onestep(const double *);
        void laplace(fftw_complex );
        void gradient(fftw_complex );
        void constant(fftw_complex ,const double *);
};

#endif //__KPSOLVER_H__
