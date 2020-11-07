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

    uid = std::uniform_int_distribution<int>(_rangel, _ranger);
    const gsl_qrng_type *qtype = nullptr;
    switch (_rtype)
    {
    case SOBOL:
        qtype = gsl_qrng_sobol;
        break;
    case HALTON:
        qtype = gsl_qrng_halton;
        break;
    case REVERSE_HALTON:
        qtype = gsl_qrng_reversehalton;
        break;
    case NIEDERREITER_2:
        qtype = gsl_qrng_niederreiter_2;
        break;
    default:
        break;
    }
    if (qtype)
    {
        qrng = gsl_qrng_alloc(qtype, 1);
        if (real_random)
        {
            int n_drop = uid(eng);
            if (n_drop < 0) n_drop = -n_drop;
            double v;
            for (int i = 0; i < n_drop; ++i)
            {
                gsl_qrng_get(qrng, &v);
            }
        }
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