#include <iostream>
#include <chrono>
#include <algorithm>
#include <list>

#include <math.h>
#include <stdlib.h>
#include <random>

#include <gsl/gsl_qrng.h>

#include "random_sampler.h"

using std::vector;
using std::list;

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

vector<unsigned> random_permutation(unsigned n)
{
    vector<unsigned> result(n);
    for (unsigned i = 0; i < n; i++)
    result[i] = i;
    std::random_shuffle(result.begin(), result.end(), randomfunc);
    return result;
}

void uniform_int_generator::init_state() 
{
    if (rangel > ranger) 
    throw std::invalid_argument("rangel is larger than ranger.");

    if (real_random) 
    {
        std::random_device rd;
        eng.seed(std::chrono::system_clock::now().time_since_epoch().count());
    }
    else 
    {
        eng.seed(0); 
    }
    if (rtype == UNIFORM) 
    {
        uid = std::uniform_int_distribution<int>(rangel, ranger);
    }
    else if (rtype == SOBOL) 
    {
        qrng = gsl_qrng_alloc(gsl_qrng_sobol, 1);
    }
    else if (rtype == HALTON)
    {
        qrng = gsl_qrng_alloc(gsl_qrng_halton, 1);
    }
    else if (rtype == REVERSE_HALTON)
    {
        qrng = gsl_qrng_alloc(gsl_qrng_reversehalton, 1);
    }
    else if (rtype == NIEDERREITER_2)
    {
        qrng = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 1);
    }
}


int uniform_int_generator::get() 
{
    switch (rtype)
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
        sample = static_cast<int>(v * (ranger - rangel) + rangel);
        break;
    case GRID:
        throw std::runtime_error("SHUFFULE not implemented.");
        break;
    }
    return sample;
}
