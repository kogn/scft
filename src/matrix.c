#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <fftw3.h>
//#include <mathlink.h>
#include <soft/utils_so3.h>
#include <math.h>

#include "wstp.h"

WSENV ep = (WSENV)0;
WSLINK lp = (WSLINK)0;


void error( WSLINK lp){
    if( WSError( lp)){
        fprintf( stderr, "Error detected by WSTP: %s.\n",
                WSErrorMessage(lp));
    }else{
        fprintf( stderr, "Error detected by this program.\n");
    }
    exit(3);
}

void deinit( void){
    if( ep) WSDeinitialize( ep);
}


void closelink( void){
    if( lp) WSClose( lp);
}

void init_and_openlink( ){
    int err;
    printf("WSINTERFACE = %d\n", WSINTERFACE);

    ep =  WSInitialize( (WSParametersPointer)0);
    if( ep == (WSENV)0) exit(1);
    atexit( deinit);

	lp = WSOpenString(ep,"-linkname \"math -mathlink\"",&err);
    //lp = WSOpenArgcArgv( ep, argc, argv, &err);
    if(lp == (WSLINK)0) exit(2);
    atexit( closelink);
}

/* Skip any packets before th   e first ReturnPacket */
/* When you can send the Wolfram System an EvaluatePacket[input], it
 *    may in general produce many packets in response, but the final
 *       packet should be ReturnPacket[output]. "Manipulating Expressions in
 *          External Programs" discusses how to handle sequences of packets and
 *             expressions whose structure you do not know in advance.  */
void skip_packets(WSLINK lp) {
    int pkt;
    while( (pkt = WSNextPacket( lp), pkt) && pkt != RETURNPKT) {
        WSNewPacket( lp);
        if (WSError( lp)) error( lp);
    }
}


void getmatrix(int bw, double alpha, double beta, double kappa, double tau,
    fftw_complex * gamma, double dt, fftw_complex *matrix, fftw_complex *matrix1)
{
  int l;
  init_and_openlink();

  WSPutFunction( lp, "EnterTextPacket", 1L);
  WSPutString(lp,"Get[\"src/matrix.m\"]");
  WSEndPacket(lp);
  for(l = 0; l<bw; l++)
  {
      WSPutFunction( lp, "EvaluatePacket", 1L);{
          WSPutFunction( lp, "matrix", 7L);{
              WSPutInteger( lp, l);
              WSPutReal64( lp, alpha);
              WSPutReal64( lp, beta);
              WSPutReal64( lp, kappa);
              WSPutReal64( lp, tau);
              WSPutReal64( lp, dt*gamma[0][0]);
              WSPutReal64( lp, dt*gamma[0][1]);
          }}WSEndPacket( lp);

    skip_packets(lp);

    double * data; 
    int *dims;
    char **heads;
    int d;
    WSGetReal64Array(lp, &data,&dims,&heads,&d);

    memcpy(matrix+totalCoeffs_so3(l),data,sizeof(fftw_complex)*(2*l+1)*(2*l+1));
    WSReleaseReal64Array(lp,data,dims,heads,d);
  }
  for(l = 0; l<bw; l++)
  {
    WSPutFunction( lp, "EvaluatePacket", 1L);
    WSPutFunction( lp, "matrix", 7L);
    WSPutInteger( lp, l);
    WSPutReal64( lp, alpha);
    WSPutReal64( lp, beta);
    WSPutReal64( lp, kappa);
    WSPutReal64( lp, tau);
    WSPutReal64( lp, dt*gamma[1][0]);
    WSPutReal64( lp, dt*gamma[1][1]);
    WSEndPacket( lp);

    skip_packets(lp);

    double * data; 
    int *dims;
    char **heads;
    int d;
    WSGetReal64Array(lp, &data,&dims,&heads,&d);

    memcpy(matrix1+totalCoeffs_so3(l),data,sizeof(fftw_complex)*(2*l+1)*(2*l+1));
    WSReleaseReal64Array(lp,data,dims,heads,d);
  }
  WSPutFunction( lp, "Exit", 0L);

  return;
}


