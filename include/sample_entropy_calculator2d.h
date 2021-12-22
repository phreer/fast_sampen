#ifndef __SAMPLE_ENTROPY_CALCULATOR2D__
#define __SAMPLE_ENTROPY_CALCULATOR2D__

#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "sample_entropy_calculator.h"
#include "global_defs.h"
#include "utils.h"

namespace sampen {
using std::vector;


template <typename T, unsigned K,
          typename =
              typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
class SampleEntropyCalculator2D : public SampleEntropyCalculator<T, K> {
  using SampleEntropyCalculator<T, K>::_data;
  using SampleEntropyCalculator<T, K>::_r;
  using SampleEntropyCalculator<T, K>::_computed;
  using SampleEntropyCalculator<T, K>::_a;
  using SampleEntropyCalculator<T, K>::_b;
  using SampleEntropyCalculator<T, K>::_output_level;
  using SampleEntropyCalculator<T, K>::_elapsed_seconds;
public:
  SampleEntropyCalculator2D(typename vector<T>::const_iterator first,
                            typename vector<T>::const_iterator last,
                            int width, int height,
                            T r,
                            OutputLevel output_level)
      : SampleEntropyCalculator<T, K>(first, last, r, output_level),
        _width(width), _height(height) {
    assert (width > K && height > K);
  }

protected:
  // Width and height of the image
  unsigned _width;
  unsigned _height;
};


template <typename T, unsigned K>
class SampleEntropyCalculator2DDirect: public SampleEntropyCalculator2D<T, K> {
  using SampleEntropyCalculator2D<T, K>::_data;
  using SampleEntropyCalculator2D<T, K>::_r;
  using SampleEntropyCalculator2D<T, K>::_computed;
  using SampleEntropyCalculator2D<T, K>::_a;
  using SampleEntropyCalculator2D<T, K>::_b;
  using SampleEntropyCalculator2D<T, K>::_output_level;
  using SampleEntropyCalculator2D<T, K>::_elapsed_seconds;
protected:
  virtual void _ComputeSampleEntropy() override;
  virtual std::string _Method() const { return "direct"; }
};


template <typename T, unsigned K>
class SampleEntropyCalculator2DSampling
    : public SampleEntropyCalculator2D<T, K> {
public:
  using SampleEntropyCalculator2D<T, K>::ComputeSampleEntropy;
  using SampleEntropyCalculator2D<T, K>::get_entropy;
  SampleEntropyCalculator2DSampling(typename vector<T>::const_iterator first,
                                    typename vector<T>::const_iterator last,
                                    int width, int height,
                                    T r, unsigned sample_size,
                                    unsigned sample_num, double real_entropy,
                                    double real_a_norm, double real_b_norm,
                                    OutputLevel output_level)
      : SampleEntropyCalculator2D<T, K>(first, last, r, output_level),
        _sample_size(sample_size), _sample_num(sample_num),
        _real_entropy(real_entropy), _real_a_norm(real_a_norm),
        _real_b_norm(real_b_norm) {}
  vector<long long> get_a_vec() {
    if (!_computed)
      ComputeSampleEntropy();
    return _a_vec;
  }
  vector<long long> get_b_vec() {
    if (!_computed)
      ComputeSampleEntropy();
    return _b_vec;
  }
  long long get_a() override {
    vector<long long> a_vec = get_a_vec();
    return std::accumulate(a_vec.cbegin(), a_vec.cend(), 0ll);
  }
  long long get_b() override {
    vector<long long> b_vec = get_b_vec();
    return std::accumulate(b_vec.cbegin(), b_vec.cend(), 0ll);
  }
  double get_a_norm() override {
    double norm =
        static_cast<double>(_sample_num * _sample_size * (_sample_size - 1));
    return get_a() / norm;
  }
  double get_b_norm() override {
    double norm =
        static_cast<double>(_sample_num * _sample_size * (_sample_size - 1));
    return get_b() / norm;
  }
  double get_err_entropy() { return get_entropy() - _real_entropy; }
  double get_err_a() { return get_a_norm() - _real_a_norm; }
  double get_err_b() { return get_b_norm() - _real_b_norm; }
  std::string get_result_str() override {
    double entropy = get_entropy();
    double error = entropy - _real_entropy;
    double rel_error = error / (entropy + 1e-8);

    std::stringstream ss;
    ss.precision(kResultDisplayPrecision);
    ss << std::scientific;
    ss << this->SampleEntropyCalculator2D<T, K>::get_result_str()
       << "\terror: " << error << "\terror (relative): " << rel_error << "\n"
       << "\terror (a): "
       << (get_a_norm() - _real_a_norm) / (_real_a_norm + 1e-8)
       << ",error (b): "
       << (get_b_norm() - _real_b_norm) / (_real_b_norm + 1e-8) << "\n";
    if (this->_output_level >= Info) {
      ss << "[INFO] sample_size: " << _sample_size
         << "\tsample_num: " << _sample_num << "\n";
      vector<long long> a_vec = get_a_vec();
      vector<long long> b_vec = get_b_vec();
      for (unsigned i = 0; i < _sample_num; ++i) {
        ss << "[INFO] "
           << "a: " << a_vec[i] << ", b: " << b_vec[i] << "\n";
      }
    }
    return ss.str();
  }

protected:
  using SampleEntropyCalculator2D<T, K>::_data;
  using SampleEntropyCalculator2D<T, K>::_r;
  using SampleEntropyCalculator2D<T, K>::_n;
  using SampleEntropyCalculator2D<T, K>::_computed;
  using SampleEntropyCalculator2D<T, K>::_a;
  using SampleEntropyCalculator2D<T, K>::_b;
  using SampleEntropyCalculator2D<T, K>::_output_level;
  using SampleEntropyCalculator2D<T, K>::_elapsed_seconds;
  const unsigned _sample_size;
  const unsigned _sample_num;
  const double _real_entropy;
  const double _real_a_norm;
  const double _real_b_norm;
  vector<long long> _a_vec, _b_vec;
};

} // namespace sampen

#endif // !__SAMPLE_ENTROPY_CALCULATOR2D__
