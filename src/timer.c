#include<sys/time.h>
#include"timer.h"
// Timer function: Measure the current time
double timer(void) {
    struct timeval Tvalue;
    struct timezone dummy;
    gettimeofday(&Tvalue, &dummy);
    double etime = (double)Tvalue.tv_sec + 1.0e-6*((double)Tvalue.tv_usec);
    return etime;
    //return omp_get_wtime();
}


