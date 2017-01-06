#ifndef __SOLVER_H__
#define __SOLVER_H__

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus
#include <fftw3.h>

#ifdef __cplusplus
}
#endif //__cplusplus
#include<string>

#ifndef NUM_THREADS
#define NUM_THREADS 8
#endif //NUM_THREADS
#ifndef DIM 
#define DIM  1
#endif

#include "transform.h"

class Solver : public Space_trans, public SO3_trans
{
    public:
        Solver(const Config &);
        ~Solver();
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
        /* friend void connect_bond_forward(TA&,TB&); */
        /* template<typename TA, typename TB> */
        /* friend void connect_bond_backward(TA&,TB&); */
        bool head_tail;

        double * phi;
        double * S[6];
        double * dist;
        double * dists;
        static double * hist_forward;
        static double * hist_backward;
    private:
        static int count;

        double kappa, tau, alpha, beta;
        fftw_complex gamma[2];

        double volume;
        double domain[DIM];

        static fftw_complex * matrix;
        static fftw_complex * matrix1;

        void tensor();

        double t, dt;

        void onestep(const double *);
        void laplace(fftw_complex);
        void gradient(fftw_complex);
        void constant(fftw_complex,const double *);

};
#endif //__SOLVER_H__

