#include <algorithm>
#include <chrono>
#include <iostream>

#include <math.h>
#include <random>
#include <stdlib.h>

#include <gsl/gsl_qrng.h>

#include "random_sampler.h"

using std::vector;

unsigned powu(unsigned x, unsigned u) {
  unsigned result = x;
  for (unsigned i = 1; i < u; i++) {
    result *= x;
  }
  return result;
}

// random generator function
int randomfunc(int j) { return rand() % j; }

void RandomIndicesSampler::init_state() {
  if (_rangel > _ranger)
    throw std::invalid_argument("_rangel is larger than _ranger.");

  if (real_random) {
    std::random_device rd;
    eng.seed(std::chrono::system_clock::now().time_since_epoch().count());
  } else {
    eng.seed(0);
  }

  uid = std::uniform_int_distribution<int>(_rangel, _ranger);
  const gsl_qrng_type *qtype = nullptr;
  switch (_rtype) {
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
  if (qtype) {
    qrng = gsl_qrng_alloc(qtype, 1);
    if (real_random) {
      _shift = uid(eng);
    } else {
      _shift = 0;
    }
  }
}

int RandomIndicesSampler::get() {
  int l = _ranger - _rangel;
  switch (_rtype) {
  case UNIFORM:
    sample = uid(eng);
    break;
  case SOBOL:
  case HALTON:
  case REVERSE_HALTON:
  case NIEDERREITER_2:
    double v;
    gsl_qrng_get(qrng, &v);
    if (real_random) {
      sample = (static_cast<int>(v * l) + _shift) % l + _rangel;
    } else {
      sample = static_cast<int>(v * l) + _rangel;
    }
    break;
  case GRID:
    MSG_ERROR(-1, "Please use class RandomIndicesSamplerSWR.\n");
  case SWR_UNIFORM:
    MSG_ERROR(-1, "Please use class RandomIndicesSamplerSWR.\n");
  }
  return sample;
}

vector<vector<unsigned>> GetSampleIndices(RandomType rtype, unsigned count,
                                          unsigned sample_size,
                                          unsigned sample_num,
                                          bool real_random) {
  if (rtype == SWR_UNIFORM || rtype == GRID) {
    return GetSampleIndicesWR(rtype, count, sample_size, sample_num,
                              real_random);
  }
  vector<vector<unsigned>> results(sample_num, vector<unsigned>(sample_size));
  RandomIndicesSampler sampler(0, count - 1, rtype, real_random);
  for (unsigned i = 0; i < results.size(); ++i) {
    for (unsigned j = 0; j < sample_size; ++j) {
      results[i][j] = static_cast<unsigned>(sampler.get());
    }
  }
  return results;
}

void RandomIndicesSamplerWR::_InitState() {
  if (_sample_size > _pop_size) {
    MSG_ERROR(-1, "The sample size (%u) is larger than population (%u).",
              _sample_size, _pop_size);
  }

  if (_real_random) {
    _eng1.seed(std::chrono::system_clock::now().time_since_epoch().count());
  } else {
    _eng1.seed(0);
  }
}

vector<vector<unsigned>> RandomIndicesSamplerWR::GetSampleArrays() {
  _uid = std::uniform_int_distribution<unsigned>();
  _urd = std::uniform_real_distribution<double>(0., 1.);
  if (_random_type == SWR_UNIFORM) {
    for (unsigned i = 0; i < _sample_num; ++i) {
      unsigned count = 0;
      unsigned count_selected = 0;
      _eng2.seed(_uid(_eng1));
      while (count < _pop_size) {
        double u = _urd(_eng2) * _pop_size / (_pop_size + 1);
        if ((_pop_size - count) * u < (_sample_size - count_selected)) {
          _samples[i][count_selected] = count;
          count_selected++;
          count++;
        } else {
          count++;
        }
      }
    }
  } else if (_random_type == GRID) {
    unsigned gap = _pop_size / _sample_size;
    unsigned offset = gap / _sample_num;
    for (unsigned i = 0; i < _sample_num; ++i) {
      unsigned index = (i * offset + _uid(_eng1)) % gap;
      for (unsigned j = 0; j < _sample_size; ++j) {
        _samples[i][j] = index;
        index += gap;
      }
    }
  }
  return _samples;
}

vector<vector<unsigned>> GetSampleIndicesWR(RandomType random_type,
                                            unsigned pop_size,
                                            unsigned sample_size,
                                            unsigned sample_num,
                                            bool real_random) {
  RandomIndicesSamplerWR sampler(pop_size, sample_size, sample_num, random_type,
                                 real_random);
  return sampler.GetSampleArrays();
}