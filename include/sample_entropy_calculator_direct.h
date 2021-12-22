#ifndef __SAMPLE_ENTROPY_CALCULATOR_DIRECT__
#define __SAMPLE_ENTROPY_CALCULATOR_DIRECT__

#include <numeric>
#include <string.h>
#include <vector>

#include "random_sampler.h"
#include "sample_entropy_calculator.h"

namespace sampen {

/**
 * Calculates an estimate of sample entropy with the raw value of r instead of
 * the variance of the estimate.
 * @param y: time series
 * @param n: the length of y
 * @param r: threshold (note that it won't be scaled with the standard deviation
 * of y)
 */
template <typename T, unsigned K>
vector<long long> _ComputeABFastDirect(const T *y, unsigned n, T r);

template <typename T, unsigned K>
vector<long long> ComputeABDirect(const vector<KDPoint<T, K + 1>> &points, T r);

template <typename T, unsigned K>
class SampleEntropyCalculatorDirect : public SampleEntropyCalculator<T, K> {
public:
  SampleEntropyCalculatorDirect(typename vector<T>::const_iterator first,
                                typename vector<T>::const_iterator last, T r,
                                OutputLevel output_level)
      : SampleEntropyCalculator<T, K>(first, last, r, output_level) {}
  std::string get_result_str() override {
    std::stringstream ss;
    ss << this->SampleEntropyCalculator<T, K>::get_result_str();
    ss << "----------------------------------------"
       << "----------------------------------------\n";
    return ss.str();
  }

protected:
  void _ComputeSampleEntropy() override;
  std::string _Method() const override { return std::string("plain direct"); }
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
class SampleEntropyCalculatorFastDirect : public SampleEntropyCalculator<T, K> {
public:
  SampleEntropyCalculatorFastDirect(typename vector<T>::const_iterator first,
                                    typename vector<T>::const_iterator last,
                                    T r, OutputLevel output_level)
      : SampleEntropyCalculator<T, K>(first, last, r, output_level) {}
  std::string get_result_str() override {
    std::stringstream ss;
    ss << this->SampleEntropyCalculator<T, K>::get_result_str();
    ss << "----------------------------------------"
       << "----------------------------------------\n";
    return ss.str();
  }

protected:
  void _ComputeSampleEntropy() override;
  std::string _Method() const override { return std::string("fast direct"); }
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
class SampleEntropyCalculatorSamplingDirect
    : public SampleEntropyCalculatorSampling<T, K> {
public:
  SampleEntropyCalculatorSamplingDirect(
      typename vector<T>::const_iterator first,
      typename vector<T>::const_iterator last, T r, unsigned sample_size,
      unsigned sample_num, double real_entropy, double real_a_norm,
      double real_b_norm, RandomType rtype, bool random_, bool presort,
      OutputLevel output_level)
      : SampleEntropyCalculatorSampling<T, K>(
            first, last, r, sample_size, sample_num, real_entropy, real_a_norm,
            real_b_norm, output_level),
        _rtype(rtype), _random(random_), _presort(presort) {}
  std::string get_result_str() override {
    std::stringstream ss;
    ss << this->SampleEntropyCalculatorSampling<T, K>::get_result_str();
    ss << "----------------------------------------"
       << "----------------------------------------\n";
    return ss.str();
  }
  double get_a_norm() override {
    return SampleEntropyCalculatorSampling<T, K>::get_a_norm() - 1. / (_n - K);
  }
  double get_b_norm() override {
    return SampleEntropyCalculatorSampling<T, K>::get_b_norm() - 1. / (_n - K);
  }

protected:
  void _ComputeSampleEntropy() override;
  std::string _Method() const override {
    std::string method_name =
        std::string("sampling direct (") + random_type_names[_rtype];
    if (_presort)
      method_name += std::string(", presort");
    method_name += std::string(")");
    return method_name;
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
  bool _presort;
};

} // namespace sampen

#endif // !__SAMPLE_ENTROPY_CALCULATOR_DIRECT__
