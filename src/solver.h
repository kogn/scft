#ifndef __SOLVER_H__
#define __SOLVER_H__

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus
#include <fftw3.h>

#ifdef __cplusplus
}
#endif //__cplusplus

#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif //NUM_THREADS
#ifndef DIM 
#define DIM  1
#endif

#include "transform.h"

class Solver : public Space_trans, public SO3_trans
{
  public:
    Solver(int, int, int[],double,double,double,double,double[]);
    Solver();
    ~Solver();

    int n_step;
    double t, dt;

    double * hist_forward;
    double * hist_backward;

    double volume;
    double domain[DIM];

    fftw_complex * matrix;
    fftw_complex * matrix1;

    fftw_complex gamma[2];

    void solve_eqn(double *);

    void init_data();
    void onestep(double *);

    double kappa, tau, alpha, beta;
    void laplace(fftw_complex);
    void gradient(fftw_complex);
    void constant(fftw_complex, double *);
};
#endif //__SOLVER_H__

