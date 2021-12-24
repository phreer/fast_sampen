#ifndef __SAMPLE_ENTROPY_CALCULATOR__
#define __SAMPLE_ENTROPY_CALCULATOR__

#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "global_defs.h"
#include "utils.h"

namespace sampen {
using std::vector;


template <typename T, unsigned K,
          typename =
              typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
class SampleEntropyCalculator {
public:
  SampleEntropyCalculator(typename vector<T>::const_iterator first,
                          typename vector<T>::const_iterator last, T r,
                          OutputLevel output_level)
      : _data(first, last), _r(r), _n(last - first),
        _output_level(output_level) {}
  virtual std::string get_result_str() {
    if (!_computed)
      ComputeSampleEntropy();
    std::stringstream ss;
    ss.precision(kResultDisplayPrecision);
    ss << "----------------------------------------"
       << "----------------------------------------\n"
       << get_method_name() << ": \n"
       << "\tsampen: " << get_entropy() << "\n"
       << "\ta (norm): " << get_a_norm() << ", b (norm): " << get_b_norm()
       << "\n"
       << "\ttime: " << std::scientific << _elapsed_seconds << "\n";
    if (_output_level >= Info) {
      std::cout << "[INFO] a: " << _a << ", b: " << _b << std::endl;
    }
    return ss.str();
  }
  double get_computation_time() {
    if (!_computed)
      this->ComputeSampleEntropy();
    return _elapsed_seconds;
  }
  double get_entropy() {
    if (!_computed)
      this->ComputeSampleEntropy();
    return ComputeSampen(static_cast<double>(get_a()),
                         static_cast<double>(get_b()), _n - K, K);
  }
  virtual long long get_a() {
    if (!_computed)
      ComputeSampleEntropy();
    return _a;
  }
  virtual long long get_b() {
    if (!_computed)
      ComputeSampleEntropy();
    return _b;
  }
  virtual double get_a_norm() {
    double norm = static_cast<double>(_n - K - 1) * (_n - K);
    return get_a() / norm;
  }
  virtual double get_b_norm() {
    double norm = static_cast<double>(_n - K - 1) * (_n - K);
    return get_b() / norm;
  }
  void ComputeSampleEntropy() {
    Timer timer;
    timer.SetStartingPointNow();
    _ComputeSampleEntropy();
    _elapsed_seconds = timer.ElapsedSeconds();
    _computed = true;
  }
  virtual std::string get_method_name() { return _Method(); }

protected:
  virtual void _ComputeSampleEntropy() = 0;
  virtual std::string _Method() const = 0;
  const vector<T> _data;
  const T _r;
  const unsigned _n;
  OutputLevel _output_level;
  long long _a, _b;
  bool _computed = false;
  double _elapsed_seconds;
};

template <typename T, unsigned K>
class SampleEntropyCalculatorSampling : public SampleEntropyCalculator<T, K> {
public:
  using SampleEntropyCalculator<T, K>::ComputeSampleEntropy;
  using SampleEntropyCalculator<T, K>::get_entropy;
  SampleEntropyCalculatorSampling(typename vector<T>::const_iterator first,
                                  typename vector<T>::const_iterator last, T r,
                                  unsigned sample_size, unsigned sample_num,
                                  double real_entropy, double real_a_norm,
                                  double real_b_norm, OutputLevel output_level)
      : SampleEntropyCalculator<T, K>(first, last, r, output_level),
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
    std::stringstream ss;
    ss.precision(kResultDisplayPrecision);
    ss << std::scientific;
    double entropy = get_entropy();
    double error = entropy - _real_entropy;
    double rel_error = error / (entropy + 1e-8);
    ss << this->SampleEntropyCalculator<T, K>::get_result_str()
       << "\terror: " << error << ", error (relative): " << rel_error << "\n"
       << "\terror (a): "
       << (get_a_norm() - _real_a_norm) / (_real_a_norm + 1e-8)
       << ", error (b): "
       << (get_b_norm() - _real_b_norm) / (_real_b_norm + 1e-8) << "\n";
    if (this->_output_level >= Info) {
      ss << "[INFO] sample_size: " << _sample_size
         << ", sample_num: " << _sample_num << "\n";
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
  using SampleEntropyCalculator<T, K>::_data;
  using SampleEntropyCalculator<T, K>::_r;
  using SampleEntropyCalculator<T, K>::_n;
  using SampleEntropyCalculator<T, K>::_computed;
  using SampleEntropyCalculator<T, K>::_a;
  using SampleEntropyCalculator<T, K>::_b;
  using SampleEntropyCalculator<T, K>::_output_level;
  using SampleEntropyCalculator<T, K>::_elapsed_seconds;
  const unsigned _sample_size;
  const unsigned _sample_num;
  const double _real_entropy;
  const double _real_a_norm;
  const double _real_b_norm;
  vector<long long> _a_vec, _b_vec;
};

} // namespace sampen

#endif // !__SAMPLE_ENTROPY_CALCULATOR__
