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
 * @param m: template length
 * of y)
 */
template <typename T>
vector<long long> _ComputeABFastDirect(const T *y, unsigned n, T r, unsigned m);

template <typename T>
vector<long long> ComputeABDirect(const vector<KDPoint<T> > &points, T r);

template <typename T>
class SampleEntropyCalculatorDirect : public SampleEntropyCalculator<T> {
public:
  std::string get_result_str() override {
    std::stringstream ss;
    ss << this->SampleEntropyCalculator<T>::get_result_str();
    ss << "----------------------------------------"
       << "----------------------------------------\n";
    return ss.str();
  }

  USING_CALCULATOR_FIELDS
protected:
  void _ComputeSampleEntropy() override;
  std::string _Method() const override { return std::string("plain direct"); }
};

template <typename T>
class SampleEntropyCalculatorFastDirect : public SampleEntropyCalculator<T> {
public:
  std::string get_result_str() override {
    std::stringstream ss;
    ss << this->SampleEntropyCalculator<T>::get_result_str();
    ss << "----------------------------------------"
       << "----------------------------------------\n";
    return ss.str();
  }

  USING_CALCULATOR_FIELDS
protected:
  void _ComputeSampleEntropy() override;
  std::string _Method() const override { return std::string("fast direct"); }
};

template <typename T>
class SampleEntropyCalculatorSamplingDirect
    : public SampleEntropyCalculatorSampling<T> {
public:
  SampleEntropyCalculatorSamplingDirect(
      typename vector<T>::const_iterator first,
      typename vector<T>::const_iterator last,
      T r, unsigned m,
      unsigned sample_size, unsigned sample_num,
      double real_entropy, double real_a_norm,
      double real_b_norm, RandomType rtype, bool random_, bool presort,
      OutputLevel output_level)
      : SampleEntropyCalculatorSampling<T>(
            first, last, r, m, sample_size, sample_num, real_entropy,
            real_a_norm, real_b_norm, output_level),
        _rtype(rtype), _random(random_), _presort(presort) {}
  SampleEntropyCalculatorSamplingDirect(
      const vector<T> &data,
      T r, unsigned m,
      unsigned sample_size, unsigned sample_num,
      double real_entropy, double real_a_norm,
      double real_b_norm, RandomType rtype, bool random_, bool presort,
      OutputLevel output_level)
      : SampleEntropyCalculatorSampling<T>(
            data.cbegin(), data.cend(), r, m, sample_size, sample_num,
            real_entropy, real_a_norm, real_b_norm, output_level),
        _rtype(rtype), _random(random_), _presort(presort) {}
  virtual std::string get_result_str() override {
    std::stringstream ss;
    ss << this->SampleEntropyCalculatorSampling<T>::get_result_str();
    ss << "----------------------------------------"
       << "----------------------------------------\n";
    return ss.str();
  }
  double get_a_norm() override {
    return SampleEntropyCalculatorSampling<T>::get_a_norm() - 1. / (_n - K);
  }
  double get_b_norm() override {
    return SampleEntropyCalculatorSampling<T>::get_b_norm() - 1. / (_n - K);
  }

  USING_SAMPLING_FIELDS
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

  RandomType _rtype;
  bool _random;
  bool _presort;
};

} // namespace sampen

#endif // !__SAMPLE_ENTROPY_CALCULATOR_DIRECT__
