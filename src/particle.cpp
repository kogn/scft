#include "particle.h"
#include "Config.h"
#include <cstdlib>
#include <cmath>
#ifndef DIM 
#define DIM 1
#endif

#ifndef NUM_THREADS
#define NUM_THREADS 8
#endif
Particle::Particle(const Config & configSettings){
    m[0] = configSettings.Read<int>("Grid_Size_x");
    domain[0] = configSettings.Read<int>("domain0");
    nB = configSettings.Read<double>("nB");
    if(DIM == 2){
        m[1] = configSettings.Read<int>("Grid_Size_y");
        domain[1] = configSettings.Read<int>("domain1");

    }

    md = 1;
    volume = 1.;
    for(int i = 0; i<DIM; i++){
        md *= m[i];
        volume *= domain[i];
    }
    phi = (double *)malloc(sizeof(double)*md);
}

Particle::~Particle(){
    free(phi);
}
void Particle::density(const double * field){
    double q = 0;
#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:q)
    for(int i = 0; i<md; i++){
        phi[i] = exp(-field[i]);
        q += phi[i];
    }
    Q = q/md;
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i = 0; i<md; i++){
        phi[i] /= Q/nB;
    }
    return;
}

void Particle::save_data(std::string filename)
{
  std::ofstream file(filename.c_str());
  for(int i = 0; i<md; i++)
  {
    file<<phi[i]<<std::endl;
  }
  file.close();
  return;
}
