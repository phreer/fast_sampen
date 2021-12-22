#include "sample_entropy_calculator_direct.h"

#include "utils.h"

namespace sampen {
template <typename T, unsigned K>
vector<long long> _ComputeABFastDirect(const T *y, unsigned n, T r) {
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


template <typename T, unsigned K>
vector<long long> ComputeABDirect(const vector<KDPoint<T, K + 1>> &points,
                                  T r) {
  const unsigned n = points.size();
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

template <typename T, unsigned K>
void SampleEntropyCalculatorDirect<T, K>::_ComputeSampleEntropy() {
  vector<KDPoint<T, K + 1>> points =
      GetKDPoints<T, K + 1>(SampleEntropyCalculator<T, K>::_data.cbegin(),
                            SampleEntropyCalculator<T, K>::_data.cend());
  vector<long long> ab =
      ComputeABDirect<T, K>(points, SampleEntropyCalculator<T, K>::_r);
  SampleEntropyCalculator<T, K>::_a = ab[0];
  SampleEntropyCalculator<T, K>::_b = ab[1];
}

template <typename T, unsigned K>
void SampleEntropyCalculatorFastDirect<T, K>::_ComputeSampleEntropy() {
  if (_n <= K) {
    std::cerr << "Data length is too short (n = " << _n;
    std::cerr << ", K = " << K << ")" << std::endl;
    exit(-1);
  }

  vector<long long> ab =
      _ComputeABFastDirect<T, K>(_data.data(), _data.size(), _r);
  _a = ab[0], _b = ab[1];
}

template <typename T, unsigned K>
void SampleEntropyCalculatorSamplingDirect<T, K>::_ComputeSampleEntropy() {
  const vector<vector<unsigned>> indices =
      GetSampleIndices(_rtype, _n - K, _sample_size, _sample_num, _random);
  _a_vec = vector<long long>(_sample_num);
  _b_vec = vector<long long>(_sample_num);
  Timer timer;
  timer.SetStartingPointNow();
  vector<vector<KDPoint<T, K + 1>>> points = GetKDPointsSample<T, K + 1>(
      _data.cbegin(), _data.cend(), indices, 1, _presort);
  timer.StopTimer();
  if (_output_level == Debug) {
    std::cout << "[INFO] Time consumed in sampling: " << timer.ElapsedSeconds()
              << " seconds. \n";
  }

  timer.SetStartingPointNow();
  for (unsigned i = 0; i < _sample_num; ++i) {
    vector<long long> ab = ComputeABDirect<T, K>(points[i], _r);
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

#define INSTANTIATE_DIRECT(K) \
  template class SampleEntropyCalculatorDirect<int, (K)>; \
  template class SampleEntropyCalculatorDirect<double, (K)>; \
  template class SampleEntropyCalculatorFastDirect<int, (K)>; \
  template class SampleEntropyCalculatorFastDirect<double, (K)>; \
  template class SampleEntropyCalculatorSamplingDirect<int, (K)>; \
  template class SampleEntropyCalculatorSamplingDirect<double, (K)>;

INSTANTIATE_DIRECT(1)
INSTANTIATE_DIRECT(2)
INSTANTIATE_DIRECT(3)
INSTANTIATE_DIRECT(4)
INSTANTIATE_DIRECT(5)
INSTANTIATE_DIRECT(6)
INSTANTIATE_DIRECT(7)
INSTANTIATE_DIRECT(8)
INSTANTIATE_DIRECT(9)
INSTANTIATE_DIRECT(10)
INSTANTIATE_DIRECT(11)

} // namespace sampen