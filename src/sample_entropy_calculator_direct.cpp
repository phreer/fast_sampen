#include "sample_entropy_calculator_direct.h"

#include "utils.h"
#include <vector>

namespace sampen {
template <typename T>
vector<long long> _ComputeABFastDirect(const T *y, unsigned n, T r, unsigned K) {
  T p[K + 1];
  long long A[K + 1];
  long long B[K + 1];
  long long *run = new long long[n];
  long long *lastrun = new long long[n];

  memset(p, 0, sizeof(p));
  memset(A, 0, sizeof(A));
  memset(B, 0, sizeof(B));
  memset(run, 0, sizeof(long long) * n);
  memset(lastrun, 0, sizeof(long long) * n);

  unsigned M1, j, nj, jj, m;
  T y1;

  unsigned M = K + 1;

  /* start running */
  for (unsigned i = 0; i < n - 1; i++) {
    nj = n - i - 1;
    y1 = y[i];
    for (jj = 0; jj < nj; jj++) {
      j = jj + i + 1;
      if (((y[j] - y1) <= r) && ((y1 - y[j]) <= r)) {
        run[jj] = lastrun[jj] + 1;
        M1 = M < run[jj] ? M : run[jj];
        for (m = 0; m < M1; m++) {
          A[m]++;
          if (j < n - 1)
            B[m]++;
        }
      } else
        run[jj] = 0;
    } /* for jj */
    for (j = 0; j < nj; j++)
      lastrun[j] = run[j];
  } /* for i */

  delete[] run;
  delete[] lastrun;

  vector<long long> result(2, 0);
  result[0] = A[M - 1];
  result[1] = B[M - 2];
  return result;
}


template <typename T>
vector<long long> ComputeABDirect(const vector<KDPoint<T> > &points,
                                  T r) {
  const unsigned n = points.size();
  if (n == 0) {
    return vector<long long>(2, 0ll);
  }
  const unsigned K = points[0].dim() - 1;
  vector<long long> results(2);
  long long a = 0LL, b = 0LL;
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = i + 1; j < n; ++j) {
      if (points[i].Within(points[j], r, K)) {
        ++b;
        T diff = points[j][K] - points[i][K];
        if (-r <= diff && diff <= r)
          ++a;
      }
    }
  }
  results[0] = a;
  results[1] = b;
  return results;
}

template <typename T>
void SampleEntropyCalculatorDirect<T>::_ComputeSampleEntropy() {
  vector<KDPoint<T> > points =
      GetKDPoints<T>(_data.cbegin(), _data.cend(), K + 1);
  vector<long long> ab = ComputeABDirect<T>(points, _r);
  _a = ab[0];
  _b = ab[1];
}

template <typename T>
void SampleEntropyCalculatorFastDirect<T>::_ComputeSampleEntropy() {
  if (_n <= K) {
    std::cerr << "Data length is too short (n = " << _n;
    std::cerr << ", K = " << K << ")" << std::endl;
    exit(-1);
  }

  vector<long long> ab =
      _ComputeABFastDirect<T>(_data.data(), _data.size(), _r, K);
  _a = ab[0], _b = ab[1];
}


template <typename T>
std::vector<long long> ComputeABSample(
    const std::vector<T> &data, const std::vector<unsigned> &sample_indices,
    const unsigned m, const T r) {
  const unsigned n0 = sample_indices.size();
  long long a = 0;
  long long b = 0;
  for (unsigned i = 0; i < n0 - 1; ++i) {
    const unsigned ii = sample_indices[i];
    for (unsigned j = i + 1; j < n0; ++j) {
      const unsigned jj = sample_indices[j];
      bool in = true;
      for (unsigned k = 0; k < m; ++k) {
        const unsigned ik = ii + k;
        const unsigned jk = jj + k;
        if (data[jk] > data[ik] + r || data[ik] > data[jk] + r) {
          in = false;
          break;
        }
      }
      if (in) {
        ++b;
        if (data[jj + m] <= data[ii + m] + r
            && data[ii + m] <= data[jj + m] + r) {
          ++a;
        }
      }
    }
  }
  std::vector<long long> result(2);
  result[0] = a;
  result[1] = b;
  return result;
}


template <typename T>
void SampleEntropyCalculatorSamplingDirect<T>::_ComputeSampleEntropy() {
  const vector<vector<unsigned> > indices =
      GetSampleIndices(_rtype, _n - K, _sample_size, _sample_num, _random);
  _a_vec = vector<long long>(_sample_num);
  _b_vec = vector<long long>(_sample_num);
  Timer timer;
  timer.SetStartingPointNow();
  timer.StopTimer();
  if (_output_level == Debug) {
    std::cout << "[INFO] Time consumed in sampling: " << timer.ElapsedSeconds()
              << " seconds.\n";
  }

  timer.SetStartingPointNow();
  for (unsigned i = 0; i < _sample_num; ++i) {
    auto ab = ComputeABSample(_data, indices[i], K, _r);
    _a_vec[i] = ab[0], _b_vec[i] = ab[1];
    _a += ab[0];
    _b += ab[1];
  }

  timer.StopTimer();
  if (_output_level == Debug) {
    std::cout << "[INFO] Time consumed in range counting: "
              << timer.ElapsedSeconds() << " seconds. \n";
  }
}

#define INSTANTIATE_DIRECT_CALCULATOR(TYPE) \
template vector<long long> _ComputeABFastDirect<TYPE>( \
    const TYPE *y, unsigned n, TYPE r, unsigned K); \
template vector<long long> ComputeABDirect<TYPE>( \
    const vector<KDPoint<TYPE> > &points, TYPE r); \
template class SampleEntropyCalculatorDirect<TYPE>; \
template class SampleEntropyCalculatorSamplingDirect<TYPE>; \
template class SampleEntropyCalculatorFastDirect<TYPE>;

INSTANTIATE_DIRECT_CALCULATOR(int)
INSTANTIATE_DIRECT_CALCULATOR(double)
} // namespace sampen