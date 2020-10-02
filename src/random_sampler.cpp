#include <iostream>
#include <chrono>
#include <algorithm>

#include <math.h>
#include <stdlib.h>
#include <random>

#include <gsl/gsl_qrng.h>

#include "random_sampler.h"

using std::vector;

unsigned powu(unsigned x, unsigned u) 
{
    unsigned result = x;
    for (unsigned i = 1; i < u; i++)
    {
        result *= x;
    }
    return result;
}

// random generator function
int randomfunc(int j)
{
    return rand() % j;
}

void RandomIndicesSampler::init_state() 
{
    if (_rangel > _ranger) 
    throw std::invalid_argument("_rangel is larger than _ranger.");

    if (real_random) 
    {
        std::random_device rd;
        eng.seed(std::chrono::system_clock::now().time_since_epoch().count());
    }
    else 
    {
        eng.seed(0); 
    }
    if (_rtype == UNIFORM) 
    {
        uid = std::uniform_int_distribution<int>(_rangel, _ranger);
    }
    else if (_rtype == SOBOL) 
    {
        qrng = gsl_qrng_alloc(gsl_qrng_sobol, 1);
    }
    else if (_rtype == HALTON)
    {
        qrng = gsl_qrng_alloc(gsl_qrng_halton, 1);
    }
    else if (_rtype == REVERSE_HALTON)
    {
        qrng = gsl_qrng_alloc(gsl_qrng_reversehalton, 1);
    }
    else if (_rtype == NIEDERREITER_2)
    {
        qrng = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 1);
    }
}


int RandomIndicesSampler::get() 
{
    switch (_rtype)
    {
    case UNIFORM:
        sample = uid(eng);
        break;
    case SOBOL:
    case HALTON:
    case REVERSE_HALTON:
    case NIEDERREITER_2:
        double v;
        gsl_qrng_get(qrng, &v);
        sample = static_cast<int>(v * (_ranger - _rangel) + _rangel);
        break;
    case GRID:
        throw std::runtime_error("SHUFFULE not implemented.");
        break;
    }
    return sample;
}

vector<vector<unsigned> > GetSampleIndices(RandomType rtype, 
                                           unsigned count, 
                                           unsigned sample_size, 
                                           unsigned sample_num, 
                                           bool random_) 
{
    vector<vector<unsigned> > results(sample_num, vector<unsigned>(sample_size)); 
    RandomIndicesSampler sampler(0, count - 1, rtype, random_); 
    for (unsigned i = 0; i < results.size(); ++i) 
    {
        for (unsigned j = 0; j < sample_size; ++j) 
        {
            results[i][j] = static_cast<unsigned>(sampler.get()); 
        }
    }
    return results; 
}