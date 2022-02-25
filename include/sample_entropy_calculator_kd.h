#ifndef __SAMPLE_ENTROPY_CALCULATOR_KD__
#define __SAMPLE_ENTROPY_CALCULATOR_KD__
#include <iostream>
#include <vector>

#include "kdtree.h"
#include "sample_entropy_calculator.h"
#include "random_sampler.h"
#include "utils.h"
#include <algorithm>

namespace sampen {
using std::vector;

// The kd tree of Mao Dong's version.
template <typename T, unsigned K> class MatchedPairsCalculatorSimpleKD {
public:
  MatchedPairsCalculatorSimpleKD(OutputLevel output_level)
      : _output_level(output_level) {}
  long long ComputeA(typename vector<T>::const_iterator first,
                     typename vector<T>::const_iterator last, T r);

private:
  OutputLevel _output_level;
};


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
class SampleEntropyCalculatorSimpleKD : public SampleEntropyCalculator<T, K> {
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

    MatchedPairsCalculatorSimpleKD<T, K> b_cal(this->_output_level);
    MatchedPairsCalculatorSimpleKD<T, K + 1> a_cal(this->_output_level);
    _b = b_cal.ComputeA(_data.cbegin(), _data.cend() - 1, _r);
    _a = a_cal.ComputeA(_data.cbegin(), _data.cend(), _r);
  }
  std::string _Method() const override { return std::string("simple kd tree"); }
  
  USING_CALCULATOR_FIELDS
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
  
  USING_CALCULATOR_FIELDS
};


// Use MatchedPairsCalculatorSampling2
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
    vector<unsigned> indices = GetSampleIndices(
        _rtype, _n - K, _sample_size, _sample_num, _random)[0];

    MatchedPairsCalculatorSampling2<T, K> b_cal(this->_output_level);
    MatchedPairsCalculatorSampling2<T, K + 1> a_cal(this->_output_level);
    _b = b_cal.ComputeA(_data.cbegin(), _data.cend() - 1, _r, indices);
    _a = a_cal.ComputeA(_data.cbegin(), _data.cend(), _r, indices);
    _a_vec = std::vector<long long>(1, _a);
    _b_vec = std::vector<long long>(1, _b);
  }
  std::string _Method() const override {
    return std::string("kd tree (Mao) sampling");
  }

  USING_SAMPLING_FIELDS
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


template <typename T, unsigned K>
class ABCalculatorSamplingLiu {
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


template <typename T, unsigned K>
class ABCalculatorRKD {
public:
  ABCalculatorRKD(OutputLevel output_level) : _output_level(output_level) {}
  vector<long long> ComputeAB(typename vector<T>::const_iterator first,
                              typename vector<T>::const_iterator last, T r);

private:
  OutputLevel _output_level;
};


template <typename T, unsigned K>
class ABCalculatorSamplingRKD {
public:
  ABCalculatorSamplingRKD(OutputLevel output_level)
      : _output_level(output_level) {}
  vector<long long> ComputeAB(typename vector<T>::const_iterator first,
                              typename vector<T>::const_iterator last,
                              T r, const vector<unsigned> &indices);

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
  
  USING_CALCULATOR_FIELDS
};


// Use MatchedPairsCalculatorSampling.
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
    vector<unsigned> indices = GetSampleIndices(
        _rtype, _n - K, _sample_size, _sample_num, _random)[0];

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

  USING_SAMPLING_FIELDS
  RandomType _rtype;
  bool _random;
};


// Use ABCalculatorSamplingLiu.
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
        _rtype(rtype), _random(random_) {
          
    if (sample_num != 1) {
      std::cerr << "Only support the parameter sample_num == 1.\n" << std::endl;
      exit(-1);
    }
  }

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
    vector<unsigned> indices = GetSampleIndices(
        _rtype, _n - K, _sample_size, _sample_num, _random)[0];

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
  
  USING_SAMPLING_FIELDS
  RandomType _rtype;
  bool _random;
};


/**
 * @brief This class calculates matched pairs like the class
 * MatchedPairCalculatorMao, except that only one kd tree is used.
 */
template <typename T, unsigned K>
class SampleEntropyCalculatorRKD : public SampleEntropyCalculator<T, K> {
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
    ABCalculatorRKD<T, K> abc(this->_output_level);
    vector<long long> result = abc.ComputeAB(_data.cbegin(), _data.cend(), _r);
    _a = result[0];
    _b = result[1];
  }
  std::string _Method() const override { return std::string("range kd tree"); }

  USING_CALCULATOR_FIELDS
};


// Use ABCalculatorSamplingRKD.
template <typename T, unsigned K>
class SampleEntropyCalculatorSamplingRKD
    : public SampleEntropyCalculatorSampling<T, K> {
public:
  SampleEntropyCalculatorSamplingRKD(typename vector<T>::const_iterator first,
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
    vector<unsigned> indices = GetSampleIndices(
        _rtype, _n - K, _sample_size, _sample_num, _random)[0];

    ABCalculatorSamplingRKD<T, K> ab_cal(_output_level);
    vector<long long> results = ab_cal.ComputeAB(_data.cbegin(), _data.cend(),
                                                 _r, indices);

    _a_vec = vector<long long>(_sample_num);
    _b_vec = vector<long long>(_sample_num);
    for (unsigned i = 0; i < _sample_num; ++i) {
      _a_vec[i] = results[2 * i];
      _b_vec[i] = results[2 * i + 1];
    }
  }
  std::string _Method() const override {
    return std::string("sampling range-kdtree (") + random_type_names[_rtype] +
           std::string(")");
  }
  USING_SAMPLING_FIELDS
  RandomType _rtype;
  bool _random;
};
} // namespace sampen

#endif // !__SAMPLE_ENTROPY_CALCULATOR_KD__
