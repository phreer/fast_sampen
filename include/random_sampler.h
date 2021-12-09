/* file: random_sampler.h
 * date: 2019-11-25
 * author: phree
 *
 * description: class used to generate random sequence
 */

#ifndef __RANDOM_SAMPLER_H__
#define __RANDOM_SAMPLER_H__

#include <random>
#include <string>
#include <vector>

#include <gsl/gsl_qrng.h>

#include "utils.h"

using std::vector;

enum RandomType {
  UNIFORM,
  SOBOL,
  HALTON,
  REVERSE_HALTON,
  NIEDERREITER_2,
  GRID,
  SWR_UNIFORM
};
static vector<std::string> random_type_names = {
    "uniform",        "sobol", "halton",     "reverse_halton",
    "NIEDERREITER_2", "GRID",  "SWR_UNIFORM"};

class RandomIndicesSampler {
public:
  /*
   * @param _rangel: minimun value
   * @param _ranger: maximum value
   */
  RandomIndicesSampler(int rangel, int ranger, RandomType rtype,
                       bool random_ = false)
      : _rangel(rangel), _ranger(ranger), _rtype(rtype), real_random(random_) {
    init_state();
  }
  ~RandomIndicesSampler() {
    if (_rtype == SOBOL || _rtype == HALTON || _rtype == REVERSE_HALTON ||
        _rtype == NIEDERREITER_2) {
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
  int _shift;
  // Whether to set seed randomly.
  bool real_random;
  void init_state();
};

vector<vector<unsigned>> GetSampleIndices(RandomType rtype, unsigned count,
                                          unsigned sample_size,
                                          unsigned sample_num,
                                          bool random_ = false);

class RandomIndicesSamplerWR {
public:
  /*
   * @param _rangel: minimun value
   * @param _ranger: maximum value
   */
  RandomIndicesSamplerWR(unsigned pop_size, unsigned sample_size,
                         unsigned sample_num, RandomType random_type,
                         bool real_random = false)
      : _pop_size(pop_size), _sample_size(sample_size), _sample_num(sample_num),
        _samples(sample_num, vector<unsigned>(sample_size)),
        _random_type(random_type), _real_random(real_random) {
    if (random_type != GRID && random_type != SWR_UNIFORM) {
      MSG_ERROR(-1, "Invalid random type for sampling without replacement.\n");
    }
    _InitState();
  }
  ~RandomIndicesSamplerWR() {}
  vector<vector<unsigned>> GetSampleArrays();

private:
  std::uniform_real_distribution<double> _urd;
  std::uniform_int_distribution<unsigned> _uid;
  std::default_random_engine _eng1;
  std::default_random_engine _eng2;
  unsigned _pop_size;
  unsigned _sample_size;
  unsigned _sample_num;
  vector<vector<unsigned>> _samples;
  // Whether to set seed randomly.
  RandomType _random_type;
  bool _real_random;
  void _InitState();
};

vector<vector<unsigned>> GetSampleIndicesWR(RandomType random_type,
                                            unsigned pop_size,
                                            unsigned sample_size,
                                            unsigned sample_num,
                                            bool real_random = false);
#endif // __RANDOM_SAMPLER_H__