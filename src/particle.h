#ifndef __PARTICLE_H__ 
#define __PARTICLE_H__ 

#ifndef DIM 
#define DIM 1 
#endif
#include "Config.h"

class Particle{
    public:
        double nB;
        int m[DIM];
        int md;
        double volume;
        double domain[DIM];
        Particle(const Config &);
        Particle();
        ~Particle();
        void save_data(std::string);
        double ptnfn();
        double Q;
        void density(const double * field);
        double * phi;
};


#endif//__PARTICLE_H__
