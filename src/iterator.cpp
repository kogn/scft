#include <cmath>
#include <iostream>
#include <fstream>

#include "iterator.hpp"
#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif //NUM_THREADS

#ifndef DIM 
#define DIM 1
#endif


Picard::Picard()
{
    alpha = 0.1;
    eps = 1e-4;
}

Anderson::Anderson()
{
    alpha = 0.1;
    eps = 1e-4;
    mk = 20;
}
SteepD::SteepD(std::string s):output_fileprefix(s)
{
  eps = 1e-4;
  steplength = 1;
}
#if DIM == 2
void SteepD::read_data2(std::string filename, double * data, int length, int times){
    std::ifstream file(filename.c_str());
    for(int i = 0; i<length/2; i++){
        file >> data[i];
        for(int j = 0; j<times; j++){
            data[i+j*length/2] = data[i];
        }
    }
    for(int i = 0; i<length/2; i++){
        file >> data[i+length*times/2];
        for(int j = 0; j<times; j++){
            data[i+j*length/2+length*times/2] = data[i+length*times/2];
        }
    }
    file.close();
    return;
}
#endif
void SteepD::read_data(std::string filename, double * data, int length){
    std::ifstream file(filename.c_str());
    for(int i = 0; i<length; i++){
        file >> data[i];
    }
    file.close();
    return;
}
void SteepD::save_data(std::string filename, double * data, int length){
    std::ofstream file(filename.c_str());
    for(int i = 0; i<length; i++){
        file << data[i] << std::endl;
    }
    file.close();
    return;
}


