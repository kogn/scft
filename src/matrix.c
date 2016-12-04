#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <fftw3.h>
#include <mathlink.h>
#include <soft/utils_so3.h>
#include <math.h>

static void init_and_openlink();



MLENV ep = (MLENV)0;
MLINK lp = (MLINK)0;

void getmatrix(int bw, double alpha, double beta, double kappa, double tau,
    fftw_complex * gamma, double dt, fftw_complex *matrix, fftw_complex *matrix1)
{
  int l;
  init_and_openlink();

  MLPutFunction( lp, "EnterTextPacket", 1L);
  MLPutString(lp,"Get[\"src/matrix.m\"]");
  MLEndPacket(lp);
  for(l = 0; l<bw; l++)
  {
    MLPutFunction( lp, "EvaluatePacket", 1L);
    MLPutFunction( lp, "matrix", 7L);
    MLPutInteger( lp, l);
    MLPutReal64( lp, alpha);
    MLPutReal64( lp, beta);
    MLPutReal64( lp, kappa);
    MLPutReal64( lp, tau);
    MLPutReal64( lp, dt*gamma[0][0]);
    MLPutReal64( lp, dt*gamma[0][1]);
    MLEndPacket( lp);

    /* skip any packets before the first ReturnPacket */
    int pkt;
    while( (pkt = MLNextPacket( lp), pkt) && pkt != RETURNPKT)
    {
      MLNewPacket( lp);
    }

    double * data; 
    int *dims;
    char **heads;
    int d;
    MLGetReal64Array(lp, &data,&dims,&heads,&d);

    memcpy(matrix+totalCoeffs_so3(l),data,sizeof(fftw_complex)*(2*l+1)*(2*l+1));
    MLReleaseReal64Array(lp,data,dims,heads,d);
  }
  for(l = 0; l<bw; l++)
  {
    MLPutFunction( lp, "EvaluatePacket", 1L);
    MLPutFunction( lp, "matrix", 7L);
    MLPutInteger( lp, l);
    MLPutReal64( lp, alpha);
    MLPutReal64( lp, beta);
    MLPutReal64( lp, kappa);
    MLPutReal64( lp, tau);
    MLPutReal64( lp, dt*gamma[1][0]);
    MLPutReal64( lp, dt*gamma[1][1]);
    MLEndPacket( lp);

    /* skip any packets before the first ReturnPacket */
    int pkt;
    while( (pkt = MLNextPacket( lp), pkt) && pkt != RETURNPKT)
    {
      MLNewPacket( lp);
    }

    double * data; 
    int *dims;
    char **heads;
    int d;
    MLGetReal64Array(lp, &data,&dims,&heads,&d);

    memcpy(matrix1+totalCoeffs_so3(l),data,sizeof(fftw_complex)*(2*l+1)*(2*l+1));
    MLReleaseReal64Array(lp,data,dims,heads,d);
  }
  MLPutFunction( lp, "Exit", 0);

  return;
}

static void deinit( void)
{
  if( ep) MLDeinitialize( ep);
}


static void closelink( void)
{
	if( lp) MLClose( lp);
}


static void init_and_openlink()
{
#if MLINTERFACE >= 3
	int err;
#else
	long err;
#endif /* MLINTERFACE >= 3 */

	ep = MLInitialize( (MLParametersPointer)0);
	if( ep == (MLENV)0) exit(1);
	atexit( deinit);
	lp = MLOpenString(ep,"-linkname \"math -mathlink\"",&err);
/*
#if MLINTERFACE < 3
	lp = MLOpenArgv( ep, argv, argv + argc, &err);
#else
	lp = MLOpenArgcArgv( ep, argc, argv, &err);
#endif
*/
	if(lp == (MLINK)0) exit(2);
	atexit( closelink);
}
