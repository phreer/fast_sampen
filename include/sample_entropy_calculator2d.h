#ifndef __SAMPLE_ENTROPY_CALCULATOR2D__
#define __SAMPLE_ENTROPY_CALCULATOR2D__

#include <array>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "random_sampler.h"
#include "utils.h"

namespace sampen {
using std::vector;

const unsigned DISPLAY_PRECISION = 6;

template <typename T,
          typename =
              typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
class SampleEntropyCalculator2D {
public:
  SampleEntropyCalculator2D(typename vector<T>::const_iterator first,
                            typename vector<T>::const_iterator last,
                            unsigned m, T r, unsigned width, unsigned height,
                            unsigned moving_step_size, unsigned dilation_factor,
                            OutputLevel output_level)
      : _data(first, last), K(m), _r(r), _width(width),
        _height(height), _moving_step_size(moving_step_size),
        _dilation_factor(dilation_factor),
        _window_size(_dilation_factor * K + 1),
        _num_steps_x((_width - _window_size + 1) / _moving_step_size),
        _num_steps_y((_width - _window_size + 1) / _moving_step_size),
        _num_templates(_num_steps_x * _num_steps_y),
        _output_level(output_level) {
    if (K == 0) {
      MSG_ERROR(-1, "Argument m must be positive.");
    }
    MSG_DEBUG("_num_steps_x: %u\n", _num_steps_x);
    MSG_DEBUG("_num_steps_y: %u\n", _num_steps_y);
  }
  virtual std::string get_result_str() {
    if (!_computed)
      ComputeSampleEntropy();
    std::stringstream ss;
    ss.precision(DISPLAY_PRECISION);
    ss << "----------------------------------------"
       << "----------------------------------------\n"
       << get_method_name() << ": \n"
       << "\tsampen2d: " << get_entropy() << "\n"
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
                         static_cast<double>(get_b()), GetNumTemplates(), K);
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
    double norm = static_cast<double>(_num_templates - 1) * (_num_templates - K);
    return get_a() / norm;
  }
  virtual double get_b_norm() {
    double norm = static_cast<double>(_num_templates - 1) * (_num_templates - K);
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
  unsigned GetNumTemplates() const {
    const unsigned window_size = _dilation_factor * K + 1;
    return (_width - window_size + 1) / _moving_step_size * \
      (_height - window_size + 1) / _moving_step_size;
  }
  virtual void _ComputeSampleEntropy() = 0;
  virtual std::string _Method() const = 0;
  bool _IsMatchedK(unsigned i, unsigned j, unsigned k, unsigned l) const {
    const unsigned end_y = i + _window_size;
    const unsigned end_x = j + _window_size;
    for (; i < end_y - 1; i += _dilation_factor, k += _dilation_factor) {
      const unsigned offset1 = i * _width;
      const unsigned offset2 = k * _width;
      for (unsigned jj = j, ll = l; jj < end_x - 1;
          jj += _dilation_factor, ll += _dilation_factor) {
        T val1 = _data[offset1 + jj];
        T val2 = _data[offset2 + ll];
        if (val1 > val2 + _r || val2 > val1 + _r) return false;
      }
    }
    return true;
  }

  bool _IsMatchedNext(unsigned i, unsigned j, unsigned k, unsigned l) const {
    const unsigned end_y = i + _window_size;
    const unsigned end_x = j + _window_size;
    for (unsigned ii = i, kk = k; ii < end_y;
        ii += _dilation_factor, kk += _dilation_factor) {
      const unsigned offset1 = ii * _width;
      const unsigned offset2 = kk * _width;
      T val1 = _data[offset1 + j + _window_size - 1];
      T val2 = _data[offset2 + l + _window_size - 1];
      if (val1 > val2 + _r || val2 > val1 + _r) return false;
    }
    const unsigned offset1 = (i + _window_size - 1) * _width;
    const unsigned offset2 = (k + _window_size - 1) * _width;
    for (; j < end_x - 1; j += _dilation_factor, l += _dilation_factor) {
      T val1 = _data[offset1 + j];
      T val2 = _data[offset2 + l];
      if (val1 > val2 + _r || val2 > val1 + _r) return false;
    }
    return true;
  }

  
  const vector<T> _data;
  const unsigned K;
  const T _r;
  const unsigned _width;
  const unsigned _height;
  const unsigned _moving_step_size;
  const unsigned _dilation_factor;

  const unsigned _window_size;
  const unsigned _num_steps_x;
  const unsigned _num_steps_y;
  const unsigned _num_templates;

  OutputLevel _output_level;

  long long _a = -1;
  long long _b = -1;
  bool _computed = false;
  double _elapsed_seconds;
};


template <typename T>
class SampleEntropyCalculator2DSampling : public SampleEntropyCalculator2D<T> {
public:
  using SampleEntropyCalculator2D<T>::ComputeSampleEntropy;
  using SampleEntropyCalculator2D<T>::_ComputeSampleEntropy;
  using SampleEntropyCalculator2D<T>::get_entropy;
  using SampleEntropyCalculator2D<T>::_data;
  using SampleEntropyCalculator2D<T>::K;
  using SampleEntropyCalculator2D<T>::_r;
  using SampleEntropyCalculator2D<T>::_width;
  using SampleEntropyCalculator2D<T>::_height;
  using SampleEntropyCalculator2D<T>::_window_size;
  using SampleEntropyCalculator2D<T>::_moving_step_size;
  using SampleEntropyCalculator2D<T>::_dilation_factor;
  using SampleEntropyCalculator2D<T>::_num_steps_x;
  using SampleEntropyCalculator2D<T>::_num_steps_y;
  using SampleEntropyCalculator2D<T>::_computed;
  using SampleEntropyCalculator2D<T>::_a;
  using SampleEntropyCalculator2D<T>::_b;
  using SampleEntropyCalculator2D<T>::_output_level;
  using SampleEntropyCalculator2D<T>::_elapsed_seconds;
  using SampleEntropyCalculator2D<T>::_num_templates;
  using SampleEntropyCalculator2D<T>::_IsMatchedK;
  using SampleEntropyCalculator2D<T>::_IsMatchedNext;
  SampleEntropyCalculator2DSampling(typename vector<T>::const_iterator first,
                                    typename vector<T>::const_iterator last,
                                    unsigned m, T r,
                                    unsigned width, unsigned height,
                                    unsigned moving_step_size,
                                    unsigned dilation_factor,
                                    unsigned sample_size, unsigned sample_num,
                                    double real_entropy, double real_a_norm,
                                    double real_b_norm, OutputLevel output_level)
      : SampleEntropyCalculator2D<T>(first, last, m, r, width, height,
                                     moving_step_size, dilation_factor,
                                     output_level),
        _sample_size(sample_size), _sample_num(sample_num),
        _real_entropy(real_entropy), _real_a_norm(real_a_norm),
        _real_b_norm(real_b_norm) {
    if (_num_templates < _sample_size) {
      MSG_ERROR(-1, "Sample size (N0) exceeds the number of templates.\n");
    }
  }
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
    ss.precision(DISPLAY_PRECISION);
    ss << std::scientific;
    double entropy = get_entropy();
    double error = entropy - _real_entropy;
    double rel_error = error / (entropy + 1e-8);
    ss << SampleEntropyCalculator2D<T>::get_result_str()
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
  const unsigned _sample_size;
  const unsigned _sample_num;
  const double _real_entropy;
  const double _real_a_norm;
  const double _real_b_norm;
  vector<long long> _a_vec, _b_vec;
};


template <typename T>
class SampleEntropyCalculator2DSamplingDirect :
    public SampleEntropyCalculator2DSampling<T> {
  using SampleEntropyCalculator2DSampling<T>::ComputeSampleEntropy;
  using SampleEntropyCalculator2DSampling<T>::get_entropy;
  using SampleEntropyCalculator2DSampling<T>::_data;
  using SampleEntropyCalculator2DSampling<T>::K;
  using SampleEntropyCalculator2DSampling<T>::_r;
  using SampleEntropyCalculator2DSampling<T>::_width;
  using SampleEntropyCalculator2DSampling<T>::_height;
  using SampleEntropyCalculator2DSampling<T>::_moving_step_size;
  using SampleEntropyCalculator2DSampling<T>::_window_size;
  using SampleEntropyCalculator2DSampling<T>::_dilation_factor;
  using SampleEntropyCalculator2DSampling<T>::_num_steps_x;
  using SampleEntropyCalculator2DSampling<T>::_num_steps_y;
  using SampleEntropyCalculator2DSampling<T>::_computed;
  using SampleEntropyCalculator2DSampling<T>::_a;
  using SampleEntropyCalculator2DSampling<T>::_b;
  using SampleEntropyCalculator2DSampling<T>::_output_level;
  using SampleEntropyCalculator2DSampling<T>::_elapsed_seconds;
  using SampleEntropyCalculator2DSampling<T>::_sample_size;
  using SampleEntropyCalculator2DSampling<T>::_sample_num;
  using SampleEntropyCalculator2DSampling<T>::_real_entropy;
  using SampleEntropyCalculator2DSampling<T>::_real_a_norm;
  using SampleEntropyCalculator2DSampling<T>::_real_b_norm;
  using SampleEntropyCalculator2DSampling<T>::_a_vec;
  using SampleEntropyCalculator2DSampling<T>::_b_vec;
  using SampleEntropyCalculator2DSampling<T>::_IsMatchedK;
  using SampleEntropyCalculator2DSampling<T>::_IsMatchedNext;
public:
  SampleEntropyCalculator2DSamplingDirect(
      typename vector<T>::const_iterator first,
      typename vector<T>::const_iterator last,
      unsigned m, T r,
      unsigned width, unsigned height,
      unsigned moving_step_size,
      unsigned dilation_factor,
      unsigned sample_size, unsigned sample_num,
      double real_entropy, double real_a_norm,
      double real_b_norm, OutputLevel output_level):
      SampleEntropyCalculator2DSampling<T>(first, last, m, r, width, height,
                                           moving_step_size, dilation_factor,
                                           sample_size, sample_num,
                                           real_entropy, real_a_norm,
                                           real_b_norm, output_level) {}
protected:
  virtual std::string _Method() const override { return "sampling direct"; }
  virtual void _ComputeSampleEntropy() override;
};

template <typename T>
class SampleEntropyCalculator2DDirect : public SampleEntropyCalculator2D<T> {
  using SampleEntropyCalculator2D<T>::_data;
  using SampleEntropyCalculator2D<T>::K;
  using SampleEntropyCalculator2D<T>::_a;
  using SampleEntropyCalculator2D<T>::_b;
  using SampleEntropyCalculator2D<T>::_r;
  using SampleEntropyCalculator2D<T>::_width;
  using SampleEntropyCalculator2D<T>::_height;;
  using SampleEntropyCalculator2D<T>::_num_steps_x;
  using SampleEntropyCalculator2D<T>::_num_steps_y;
  using SampleEntropyCalculator2D<T>::_moving_step_size;
  using SampleEntropyCalculator2D<T>::_dilation_factor;;
  using SampleEntropyCalculator2D<T>::_window_size;;
  using SampleEntropyCalculator2D<T>::_IsMatchedK;
  using SampleEntropyCalculator2D<T>::_IsMatchedNext;
public:
  SampleEntropyCalculator2DDirect(typename vector<T>::const_iterator first,
                            typename vector<T>::const_iterator last,
                            unsigned m, T r, unsigned width, unsigned height,
                            unsigned moving_step_size, unsigned dilation_factor,
                            OutputLevel output_level)
      : SampleEntropyCalculator2D<T>(first, last, m, r, width, height,
                                     moving_step_size, dilation_factor,
                                     output_level) {}
protected:
  virtual std::string _Method() const override { return "direct"; }
  virtual void _ComputeSampleEntropy() override;
};


// Implementation
template <typename T>
void SampleEntropyCalculator2DDirect<T>::_ComputeSampleEntropy() {
  long long a = 0;
  long long b = 0;
  const unsigned end_x = _width - _window_size + 1;
  const unsigned end_y = _height - _window_size + 1;
  // Compare template (i, j) and template (k, l).
  for (unsigned i = 0; i < end_y; i += _moving_step_size) {
    for (unsigned j = 0; j < end_x; j += _moving_step_size) {
      for (unsigned k = i + _moving_step_size; k < end_y;
          k += _moving_step_size) {
        if (_IsMatchedK(i, j, k, j)) {
          b++;
          if (_IsMatchedNext(i, j, k, j)) {
            a++;
          }
        }
      }
      for (unsigned k = 0; k < end_y; k += _moving_step_size) {
        for (unsigned l = j + _moving_step_size; l < end_x;
            l += _moving_step_size) {
          if (_IsMatchedK(i, j, k, l)) {
            b++;
            if (_IsMatchedNext(i, j, k, l)) {
              a++;
            }
          }        
        }
      }
    }
  }
  _a = a;
  _b = b;
}


template <typename T>
void SampleEntropyCalculator2DSamplingDirect<T>::_ComputeSampleEntropy() {
  unsigned num_templates = _num_steps_x * _num_steps_y;
  RandomIndicesSamplerWR sampler(num_templates, _sample_size, _sample_num,
                                 RandomType::SWR_UNIFORM);
  const auto random_indices_array = sampler.GetSampleArrays();
  _a_vec.clear();
  _b_vec.clear();
  for (unsigned sample_i = 0; sample_i < random_indices_array.size(); ++sample_i) {
    const auto& random_indices = random_indices_array[sample_i];
    std::vector<std::array<unsigned, 2> > random_indices_xy;
    for (auto index: random_indices) {
      random_indices_xy.push_back({_moving_step_size * (index / _num_steps_x),
                                   _moving_step_size * (index % _num_steps_x)});
    }
    long long a = 0;
    long long b = 0;
    for (unsigned i = 0; i < _sample_size - 1; ++i) {
      const auto& index_xy = random_indices_xy[i];
      for (unsigned j = i + 1; j < _sample_size; ++j) {
        const auto& index_xy2 = random_indices_xy[j];
        if (_IsMatchedK(index_xy[0], index_xy[1], index_xy2[0], index_xy2[1])) {
          ++b;
          if (_IsMatchedNext(index_xy[0], index_xy[1], index_xy2[0], index_xy2[1])) {
            ++a;
          }
        }
      }
    }
    _a_vec.push_back(a);
    _b_vec.push_back(b);
  }
}
} // namespace sampen

#endif // !__SAMPLE_ENTROPY_CALCULATOR2D__
