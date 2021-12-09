#ifndef __SAMPLE_ENTROPY_CALCULATOR_KD__
#define __SAMPLE_ENTROPY_CALCULATOR_KD__
#include <iostream>
#include <vector>

#include "kdtree.h"
#include "random_sampler.h"
#include "sample_entropy_calculator.h"
#include "utils.h"
#include <algorithm>

namespace sampen {
using std::vector;

// The kd tree of Mao Dong's version.
template <typename T, unsigned K> class MatchedPairsCalculatorMao {
public:
  MatchedPairsCalculatorMao(OutputLevel output_level)
      : _output_level(output_level) {}
  long long ComputeA(typename vector<T>::const_iterator first,
                     typename vector<T>::const_iterator last, T r);

private:
  OutputLevel _output_level;
};

template <typename T, unsigned K> class MatchedPairsCalculatorSampling {
public:
  MatchedPairsCalculatorSampling(OutputLevel output_level)
      : _output_level(output_level) {}
  vector<long long> ComputeA(typename vector<T>::const_iterator first,
                             typename vector<T>::const_iterator last,
                             unsigned sample_num,
                             const vector<unsigned> &indices, T r);

private:
  OutputLevel _output_level;
};

template <typename T, unsigned K> class MatchedPairsCalculatorSampling2 {
public:
  MatchedPairsCalculatorSampling2(OutputLevel output_level)
      : _output_level(output_level) {}
  long long ComputeA(typename vector<T>::const_iterator first,
                     typename vector<T>::const_iterator last, T r,
                     vector<unsigned> &indices);

private:
  OutputLevel _output_level;
};

template <typename T, unsigned K>
class SampleEntropyCalculatorMao : public SampleEntropyCalculator<T, K> {
public:
  using SampleEntropyCalculator<T, K>::SampleEntropyCalculator;
  std::string get_result_str() override {
    std::stringstream ss;
    ss << this->SampleEntropyCalculator<T, K>::get_result_str();
    ss << "----------------------------------------"
       << "----------------------------------------\n";
    return ss.str();
  }

protected:
  void _ComputeSampleEntropy() override {
    if (_n <= K) {
      std::cerr << "Data length is too short (n = " << _n;
      std::cerr << ", K = " << K << ")" << std::endl;
      exit(-1);
    }

    MatchedPairsCalculatorMao<T, K> b_cal(this->_output_level);
    MatchedPairsCalculatorMao<T, K + 1> a_cal(this->_output_level);
    _b = b_cal.ComputeA(_data.cbegin(), _data.cend() - 1, _r);
    _a = a_cal.ComputeA(_data.cbegin(), _data.cend(), _r);
  }
  std::string _Method() const override { return std::string("kd tree (Mao)"); }
  using SampleEntropyCalculator<T, K>::_data;
  using SampleEntropyCalculator<T, K>::_r;
  using SampleEntropyCalculator<T, K>::_n;
  using SampleEntropyCalculator<T, K>::_computed;
  using SampleEntropyCalculator<T, K>::_a;
  using SampleEntropyCalculator<T, K>::_b;
  using SampleEntropyCalculator<T, K>::_output_level;
  using SampleEntropyCalculator<T, K>::_elapsed_seconds;
};

template <typename T, unsigned K>
class SampleEntropyCalculatorSamplingMao
    : public SampleEntropyCalculatorSampling<T, K> {
public:
  using SampleEntropyCalculatorSampling<T, K>::SampleEntropyCalculatorSampling;
  SampleEntropyCalculatorSamplingMao(typename vector<T>::const_iterator first,
                                     typename vector<T>::const_iterator last,
                                     T r, unsigned sample_size,
                                     unsigned sample_num, double real_entropy,
                                     double real_a_norm, double real_b_norm,
                                     RandomType rtype, bool random_,
                                     OutputLevel output_level)
      : SampleEntropyCalculatorSampling<T, K>(
            first, last, r, sample_size, sample_num, real_entropy, real_a_norm,
            real_b_norm, output_level),
        _rtype(rtype), _random(random_) {
    if (sample_num != 1) {
      std::cerr << "Only support the parameter sample_num == 1.\n" << std::endl;
      exit(-1);
    }
  }
  std::string get_result_str() override {
    std::stringstream ss;
    ss << this->SampleEntropyCalculatorSampling<T, K>::get_result_str();
    ss << "----------------------------------------"
       << "----------------------------------------\n";
    return ss.str();
  }

protected:
  void _ComputeSampleEntropy() override {
    if (_n <= K) {
      std::cerr << "Data length is too short (n = " << _n;
      std::cerr << ", K = " << K << ")" << std::endl;
      exit(-1);
    }
    RandomIndicesSamplerWR sampler(_n - K, _sample_size, 1, _rtype, _random);
    std::vector<unsigned> sample_indices = sampler.GetSampleArrays()[0];
    MatchedPairsCalculatorSampling2<T, K> b_cal(this->_output_level);
    MatchedPairsCalculatorSampling2<T, K + 1> a_cal(this->_output_level);
    _b = b_cal.ComputeA(_data.cbegin(), _data.cend() - 1, _r, sample_indices);
    _a = a_cal.ComputeA(_data.cbegin(), _data.cend(), _r, sample_indices);
    _a_vec = std::vector<long long>(1, _a);
    _b_vec = std::vector<long long>(1, _b);
  }
  std::string _Method() const override {
    return std::string("kd tree (Mao) sampling");
  }
  using SampleEntropyCalculatorSampling<T, K>::_data;
  using SampleEntropyCalculatorSampling<T, K>::_r;
  using SampleEntropyCalculatorSampling<T, K>::_n;
  using SampleEntropyCalculatorSampling<T, K>::_computed;
  using SampleEntropyCalculatorSampling<T, K>::_a;
  using SampleEntropyCalculatorSampling<T, K>::_b;
  using SampleEntropyCalculatorSampling<T, K>::_output_level;
  using SampleEntropyCalculatorSampling<T, K>::_elapsed_seconds;
  using SampleEntropyCalculatorSampling<T, K>::_sample_size;
  using SampleEntropyCalculatorSampling<T, K>::_sample_num;
  using SampleEntropyCalculatorSampling<T, K>::_a_vec;
  using SampleEntropyCalculatorSampling<T, K>::_b_vec;
  RandomType _rtype;
  bool _random;
};

template <typename T, unsigned K> class ABCalculatorLiu {
public:
  ABCalculatorLiu(OutputLevel output_level) : _output_level(output_level) {}
  vector<long long> ComputeAB(typename vector<T>::const_iterator first,
                              typename vector<T>::const_iterator last, T r);

private:
  OutputLevel _output_level;
};

template <typename T, unsigned K> class ABCalculatorSamplingLiu {
public:
  ABCalculatorSamplingLiu(OutputLevel output_level)
      : _output_level(output_level) {}
  vector<long long> ComputeAB(typename vector<T>::const_iterator first,
                              typename vector<T>::const_iterator last,
                              unsigned sample_num,
                              const vector<unsigned> &indices, T r);

private:
  OutputLevel _output_level;
};

/**
 * @brief This class calculates matched pairs like the class
 * MatchedPairCalculatorMao, except that only one kd tree is used.
 */
template <typename T, unsigned K>
class SampleEntropyCalculatorLiu : public SampleEntropyCalculator<T, K> {
public:
  using SampleEntropyCalculator<T, K>::SampleEntropyCalculator;
  std::string get_result_str() override {
    std::stringstream ss;
    ss << this->SampleEntropyCalculator<T, K>::get_result_str();
    ss << "----------------------------------------"
       << "----------------------------------------\n";
    return ss.str();
  }

protected:
  void _ComputeSampleEntropy() override {
    if (_n <= K) {
      std::cerr << "Data length is too short (n = " << _n;
      std::cerr << ", K = " << K << ")" << std::endl;
      exit(-1);
    }
    ABCalculatorLiu<T, K> abc(this->_output_level);
    vector<long long> result = abc.ComputeAB(_data.cbegin(), _data.cend(), _r);
    _a = result[0];
    _b = result[1];
  }
  std::string _Method() const override { return std::string("kd tree (Liu)"); }

  using SampleEntropyCalculator<T, K>::_data;
  using SampleEntropyCalculator<T, K>::_r;
  using SampleEntropyCalculator<T, K>::_n;
  using SampleEntropyCalculator<T, K>::_computed;
  using SampleEntropyCalculator<T, K>::_a;
  using SampleEntropyCalculator<T, K>::_b;
  using SampleEntropyCalculator<T, K>::_output_level;
  using SampleEntropyCalculator<T, K>::_elapsed_seconds;
};

template <typename T, unsigned K>
class SampleEntropyCalculatorSamplingKDTree
    : public SampleEntropyCalculatorSampling<T, K> {
public:
  SampleEntropyCalculatorSamplingKDTree(
      typename vector<T>::const_iterator first,
      typename vector<T>::const_iterator last, T r, unsigned sample_size,
      unsigned sample_num, double real_entropy, double real_a_norm,
      double real_b_norm, RandomType rtype, bool random_,
      OutputLevel output_level)
      : SampleEntropyCalculatorSampling<T, K>(
            first, last, r, sample_size, sample_num, real_entropy, real_a_norm,
            real_b_norm, output_level),
        _rtype(rtype), _random(random_) {}
  std::string get_result_str() override {
    std::stringstream ss;
    ss << this->SampleEntropyCalculatorSampling<T, K>::get_result_str();
    if (_output_level) {
      ss.precision(4);
      ss << std::scientific;
      ss << "[INFO] "
         << "random_type: " << random_type_names[_rtype] << "\n"
         << "[INFO] "
         << "random: " << _random << "\n";
    }
    ss << "----------------------------------------"
       << "----------------------------------------\n";
    return ss.str();
  }

protected:
  void _ComputeSampleEntropy() override {
    if (_n <= K) {
      std::cerr << "Data length is too short (n = " << _n;
      std::cerr << ", K = " << K << ")" << std::endl;
      exit(-1);
    }
    RandomIndicesSampler uig(0, _n - 2, _rtype, _random);
    vector<unsigned> indices(_sample_size * _sample_num);
    for (unsigned i = 0; i < _sample_num; ++i) {
      for (unsigned j = 0; j < _sample_size; ++j) {
        indices[i * _sample_size + j] = uig.get();
      }
    }

    MatchedPairsCalculatorSampling<T, K> b_cal(this->_output_level);
    MatchedPairsCalculatorSampling<T, K + 1> a_cal(this->_output_level);
    _a_vec =
        a_cal.ComputeA(_data.cbegin(), _data.cend(), _sample_num, indices, _r);
    _b_vec =
        b_cal.ComputeA(_data.cbegin(), _data.cend(), _sample_num, indices, _r);
  }

  std::string _Method() const override {
    return std::string("sampling kd tree (Mao, ") + random_type_names[_rtype] +
           std::string(")");
  }

  using SampleEntropyCalculatorSampling<T, K>::_data;
  using SampleEntropyCalculatorSampling<T, K>::_r;
  using SampleEntropyCalculatorSampling<T, K>::_n;
  using SampleEntropyCalculatorSampling<T, K>::_computed;
  using SampleEntropyCalculatorSampling<T, K>::_a;
  using SampleEntropyCalculatorSampling<T, K>::_b;
  using SampleEntropyCalculatorSampling<T, K>::_output_level;
  using SampleEntropyCalculatorSampling<T, K>::_elapsed_seconds;
  using SampleEntropyCalculatorSampling<T, K>::_sample_size;
  using SampleEntropyCalculatorSampling<T, K>::_sample_num;
  using SampleEntropyCalculatorSampling<T, K>::_a_vec;
  using SampleEntropyCalculatorSampling<T, K>::_b_vec;
  RandomType _rtype;
  bool _random;
};

template <typename T, unsigned K>
class SampleEntropyCalculatorSamplingLiu
    : public SampleEntropyCalculatorSampling<T, K> {
public:
  SampleEntropyCalculatorSamplingLiu(typename vector<T>::const_iterator first,
                                     typename vector<T>::const_iterator last,
                                     T r, unsigned sample_size,
                                     unsigned sample_num, double real_entropy,
                                     double real_a_norm, double real_b_norm,
                                     RandomType rtype, bool random_,
                                     OutputLevel output_level)
      : SampleEntropyCalculatorSampling<T, K>(
            first, last, r, sample_size, sample_num, real_entropy, real_a_norm,
            real_b_norm, output_level),
        _rtype(rtype), _random(random_) {}
  std::string get_result_str() override {
    std::stringstream ss;
    ss << this->SampleEntropyCalculatorSampling<T, K>::get_result_str();
    if (_output_level >= Info) {
      ss.precision(4);
      ss << std::scientific;
      ss << "[INFO] "
         << "random_type: " << random_type_names[_rtype] << "\n"
         << "[INFO] "
         << "random: " << _random << "\n";
    }
    ss << "----------------------------------------"
       << "----------------------------------------\n";
    return ss.str();
  }

protected:
  void _ComputeSampleEntropy() override {
    if (_n <= K) {
      std::cerr << "Data length is too short (n = " << _n;
      std::cerr << ", K = " << K << ")" << std::endl;
      exit(-1);
    }
    RandomIndicesSampler uig(0, _n - 2, _rtype, _random);
    vector<unsigned> indices(_sample_size * _sample_num);
    for (unsigned i = 0; i < _sample_num; ++i) {
      for (unsigned j = 0; j < _sample_size; ++j) {
        indices[i * _sample_size + j] = uig.get();
      }
    }

    ABCalculatorSamplingLiu<T, K> ab_cal(_output_level);
    vector<long long> results = ab_cal.ComputeAB(_data.cbegin(), _data.cend(),
                                                 _sample_num, indices, _r);

    _a_vec = vector<long long>(_sample_num);
    _b_vec = vector<long long>(_sample_num);
    for (unsigned i = 0; i < _sample_num; ++i) {
      _a_vec[i] = results[2 * i];
      _b_vec[i] = results[2 * i + 1];
    }
  }
  std::string _Method() const override {
    return std::string("sampling kd tree (liu, ") + random_type_names[_rtype] +
           std::string(")");
  }

  using SampleEntropyCalculatorSampling<T, K>::_data;
  using SampleEntropyCalculatorSampling<T, K>::_r;
  using SampleEntropyCalculatorSampling<T, K>::_n;
  using SampleEntropyCalculatorSampling<T, K>::_computed;
  using SampleEntropyCalculatorSampling<T, K>::_a;
  using SampleEntropyCalculatorSampling<T, K>::_b;
  using SampleEntropyCalculatorSampling<T, K>::_output_level;
  using SampleEntropyCalculatorSampling<T, K>::_elapsed_seconds;
  using SampleEntropyCalculatorSampling<T, K>::_sample_size;
  using SampleEntropyCalculatorSampling<T, K>::_sample_num;
  using SampleEntropyCalculatorSampling<T, K>::_a_vec;
  using SampleEntropyCalculatorSampling<T, K>::_b_vec;
  RandomType _rtype;
  bool _random;
};

} // namespace sampen

#endif // !__SAMPLE_ENTROPY_CALCULATOR_KD__
