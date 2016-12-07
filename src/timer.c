#include<sys/time.h>
#include"timer.h"
static int start_time= 0;
// Timer function: Measure the current time
double timer(void) {
    struct timeval Tvalue;
    struct timezone dummy;
    gettimeofday(&Tvalue, &dummy);
    double etime = (double)Tvalue.tv_sec + 1.0e-6*((double)Tvalue.tv_usec);
    if(start_time == 0) start_time = etime;
    return etime - start_time;
    //return omp_get_wtime();
}

