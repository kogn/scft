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
SteepD::SteepD(std::string s):output_filedir(s)
{
  eps = 1e-4;
  steplength = 1;
}
