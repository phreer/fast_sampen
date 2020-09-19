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

enum RandomType {UNIFORM, SOBOL, HALTON, REVERSE_HALTON, NIEDERREITER_2, GRID};

class RandomIndicesSampler
{
public:
    /*
     * @param _rangel: minimun value
     * @param _ranger: maximum value 
     */
    RandomIndicesSampler(
        int rangel, int ranger, RandomType rtype, bool random_ = false):
            _rangel(rangel), _ranger(ranger), _rtype(rtype), 
            real_random(random_)
    {
        init_state();
    }
    ~RandomIndicesSampler() 
    {
        if (_rtype == SOBOL || _rtype == HALTON || _rtype == REVERSE_HALTON 
            || _rtype == NIEDERREITER_2) 
        {
            gsl_qrng_free(qrng); 
        }
    }
    int get();
private:
    int _rangel, _ranger;
    RandomType _rtype;
    std::uniform_int_distribution<int> uid;
    std::default_random_engine eng;
    gsl_qrng *qrng;
    int sample;
    // Whether to set seed randomly.
    bool real_random;
    void init_state();
};

vector<vector<unsigned> > GetSampleIndices(RandomType rtype, 
                                           unsigned count, 
                                           unsigned sample_size, 
                                           unsigned sample_num, 
                                           bool random_ = false); 
#endif // __RANDOM_SAMPLER_H__