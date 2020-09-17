/* file: random_sampler.h
 * date: 2019-11-25
 * author: phree
 *
 * description: class used to generate random sequence
 */

#ifndef __RANDOM_SAMPLER_H__
#define __RANDOM_SAMPLER_H__

#include <vector>
#include <random>
#include <gsl/gsl_qrng.h>

#include "utils.h"

using std::vector;

class uniform_int_generator
{
public:
    enum random_type {UNIFORM, SOBOL, HALTON, REVERSE_HALTON, NIEDERREITER_2, GRID};
    /*
     * @param _rangel: minimun value
     * @param _ranger: maximum value 
     */
    uniform_int_generator(
        int _rangel, int _ranger, random_type _rtype, bool _random = false): 
            rangel(_rangel), ranger(_ranger), rtype(_rtype), 
            real_random(_random)
    {
        init_state();
    }
    ~uniform_int_generator() 
    {
        if (rtype == SOBOL || rtype == HALTON || rtype == REVERSE_HALTON 
            || rtype == NIEDERREITER_2) 
        {
            gsl_qrng_free(qrng); 
        }
    }
    int get();
private:
    int rangel, ranger;
    const uniform_int_generator::random_type rtype;
    std::uniform_int_distribution<int> uid;
    std::default_random_engine eng;
    gsl_qrng *qrng;
    int sample;
    // Whether to set seed randomly.
    bool real_random;
    void init_state();
};

#endif // __RANDOM_SAMPLER_H__
