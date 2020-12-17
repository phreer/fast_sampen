#include <iostream>
#include <chrono>
#include <algorithm>

#include <math.h>
#include <stdlib.h>
#include <random>

#include <gsl/gsl_qrng.h>

#include "random_sampler.h"
#include <time.h>
#include <cstdlib>

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
            _shift = uid(eng);
        }
        else
        {
            _shift = 0;
        }
    }
}


int RandomIndicesSampler::get() 
{
    int l = _ranger - _rangel;
    switch (_rtype)
    {
    case UNIFORM:
        //sample = uid(eng);
        break;
    case SOBOL:
    case HALTON:
    case REVERSE_HALTON:
    case NIEDERREITER_2:
        double v;
        gsl_qrng_get(qrng, &v);
        if (real_random)
        {
            sample = (static_cast<int>(v * l) + _shift) % l + _rangel;
        }
        else
        {
            sample = static_cast<int>(v * l) + _rangel;
        }
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

//以下是修改的
void RandomIndicesSamplerWR::init_state() 
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
    case SWR:
        N = _ranger;
        n = _sample_quantity;
        t = 0;
        m = 0;
        //qtype = gsl_qrng_swr;
        break;
    default:
        break;
    }
    if (qtype)
    {
        qrng = gsl_qrng_alloc(qtype, 1);
        if (real_random)
        {
            _shift = uid(eng);
        }
        else
        {
            _shift = 0;
        }
    }
}


vector<int> RandomIndicesSamplerWR::get() 
{
    int l = _ranger - _rangel;
    switch (_rtype)
    {
    case SWR:
        srand(time(0));
        while(t <= N)
        {
            double U = rand()%100/(double)101;
            if((N - t)*U < (n - m))
            {
                sample.push_back(t);
                m++;
                t++;
            }
            else
            {
                t++;
            }
        }
        break;
    }
    return sample;
}

vector<unsigned> GetSampleIndicesWR(RandomType rtype, 
                                           unsigned count, 
                                           unsigned sample_size, 
                                           unsigned sample_num, 
                                           bool random_) 
{
    //vector<vector<unsigned> > results(sample_num, vector<unsigned>(sample_size)); 
    RandomIndicesSamplerWR sampler(0, count - 1, sample_size*sample_num, rtype, random_);
    vector<unsigned> results(sample_size*sample_num);
    vector<int> temp = sampler.get();
    int num = 0;
    for (unsigned i = 0; i < sample_size*sample_num; ++i) 
    {
        results[i] = static_cast<unsigned>(temp[num++]); 
    }
    return results; 
}