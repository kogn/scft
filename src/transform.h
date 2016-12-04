#ifndef __TRANSFORM_H__
#define __TRANSFORM_H__

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

#include <fftw3.h>

#ifdef __cplusplus
}
#endif //__cplusplus

class Data 
{
  public:
    int bw, n, n3, m, md, n_coeff;
    fftw_complex * realdata;
    fftw_complex * spectdata;
    Data(int, int);
    virtual ~Data();
};

class Space_trans:virtual public Data
{
  public:
    fftw_plan pf, pi;
    void for_space();
    void inv_space();

    Space_trans(int, int);
    ~Space_trans();
};

class SO3_trans : virtual public Data
{
  public:
    void for_so3();
    void inv_so3();
    SO3_trans(int, int);
    ~SO3_trans();
    double * weights;
  private:
    fftw_plan p1, p2;
    fftw_complex ** workspace_cx;
    fftw_complex ** workspace_cx2;
    double ** workspace_re;
    int wignerSpace;
    double * wigners;
    double * wignersTrans;
};

#endif //__TRANSFORM_H__
