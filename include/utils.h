/**
 * @file utils.h
 * @author Liu Weifeng
 * @date 2020/05/10
 *
 * @brief Some utility functions and class.
 *
 * @details
 */
#ifndef __UTILS__
#define __UTILS__
#include <algorithm>
#include <assert.h>
#include <chrono>
#include <fstream>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#ifdef DEBUG
#include <iostream>
#endif

#include "kdpoint.h"
#include "time.h"

#define MSG(fd, prefix, fmt, ...)                                              \
  fprintf(fd, "[%s] ", (prefix));                                              \
  fprintf(fd, (fmt), ##__VA_ARGS__);

#define MSG_ERROR(error_code, fmt, ...)                                        \
  MSG(stderr, "ERROR", (fmt), ##__VA_ARGS__);                                  \
  MSG(stderr, "ERROR", "File: %s, line: %d\n", __FILE__, __LINE__);            \
  exit((error_code));

#ifdef DEBUG
#define MSG_DEBUG(fmt, ...) MSG(stdout, "DEBUG", (fmt), ##__VA_ARGS__);
#else
#define MSG_DEBUG(fmt, ...)                                                    \
  {}
#endif

#define MSG_INFO(fmt, ...) MSG(stdout, "INFO", (fmt), ##__VA_ARGS__);

namespace sampen {
using std::cerr;
using std::ifstream;
using std::string;
using std::vector;

enum OutputLevel {
  // Only show essential result, i.e., necessary results.
  Silent,
  // Show helpful informations, illustrating the status of functions.
  Info,
  // Show anything that may be useful for debugging.
  Debug
};

struct Bounds {
  Bounds(size_t n) : lower_bounds(n), upper_bounds(n) {}
  vector<unsigned> lower_bounds;
  vector<unsigned> upper_bounds;
};

template <typename T, unsigned K> struct Range {
  T lower_ranges[K];
  T upper_ranges[K];
};

vector<unsigned> GetInverseMap(const vector<unsigned> &map);

double ComputeSampen(double A, double B, unsigned N, unsigned m);

/**
 * @brief Read data from file.
 *
 * @param filename: The name of the file to read.
 * @return Read data.
 */
template <typename T>
vector<T> ReadData(std::string filename, std::string input_type = "simple",
                   unsigned n = 0);

template <typename T> double ComputeVariance(const vector<T> &data);

template <typename T> T ComputeSum(const vector<T> &data);

/**
 * @brief Merge repeated points by setting count.
 *
 * @param[out] points: The k-dimensional points that has been sorted.
 * @note Auxiliary points will be disabled.
 */

template <typename T, unsigned K>
void MergeRepeatedPoints(vector<KDPoint<T, K>> &points,
                         const vector<unsigned> &rank2index);

template <typename T, unsigned K>
vector<KDPoint<T, K>> GetKDPoints(typename vector<T>::const_iterator first,
                                  typename vector<T>::const_iterator last,
                                  int count = 1);

template <typename T, unsigned K>
vector<vector<KDPoint<T, K>>>
GetKDPointsSample(typename vector<T>::const_iterator first,
                  typename vector<T>::const_iterator last,
                  const vector<vector<unsigned>> &indices, int count,
                  bool presort = false, OutputLevel output_level = Silent);

template <typename T, unsigned K>
void MergeRepeatedPoints(vector<KDPoint<T, K>> &points,
                         const vector<unsigned> &rank2index);

template <typename T, unsigned K>
void CloseAuxiliaryPoints(vector<KDPoint<T, K>> &points,
                          const vector<unsigned> &rank2index);

/*
 * @brief Maps the points to grids.
 *
 * For each point p, p[i] is mapped to the rank of p[i] for i = 1, 2, ..., K-1
 */
template <typename T, unsigned K>
vector<KDPoint<unsigned, K - 1>> Map2Grid(const vector<KDPoint<T, K>> &points,
                                          const vector<unsigned> &rank2index,
                                          bool skip_nocount = true);

/**
 * @brief Get the bounds of indices such that within the bound the values
 * are within the threshold r.
 *
 * @param data: The data to find lower bounds and upper bounds.
 * @note The data must have been sorted in ascending lexicographic order.
 * @param r: The threshold.
 * @return The bounds of the indices.
 */
template <typename T, unsigned K>
Bounds GetRankBounds(const vector<KDPoint<T, K>> &points, T r);

/*
 * @brief Given a point (in grid), get the bound.
 */
template <unsigned K>
Range<unsigned, K> GetHyperCube(const KDPoint<unsigned, K> &point,
                                const Bounds &bounds);

class ArgumentParser {
public:
  ArgumentParser(int argc, char *argv[]) : arg_list(argv, argv + argc) {}
  string getArg(const string &arg);
  bool isOption(const string &opt);
  long getArgLong(const string &arg, long default_);
  double getArgDouble(const string &arg, double default_);
  std::vector<int> getArgIntArray(const string &arg);

private:
  vector<string> arg_list;
};

/**
 * Timer class for evaluating the time elapsed from a starting point.
 */
class SysTimer {
public:
  SysTimer() { SetStartingPointNow(); }
  void SetStartingPointNow() {
    _starting_point = std::chrono::system_clock::now();
    _runing = true;
  }
  void StopTimer() {
    if (_runing) {
      _end_point = std::chrono::system_clock::now();
      _runing = false;
    }
  }
  double ElapsedMilliseconds() {
    std::chrono::time_point<std::chrono::system_clock> end_point;

    if (_runing) {
      end_point = std::chrono::system_clock::now();
    } else {
      end_point = _end_point;
    }

    return std::chrono::duration_cast<std::chrono::milliseconds>(
               end_point - _starting_point)
        .count();
  }
  double ElapsedSeconds() { return ElapsedMilliseconds() / 1000.; }

private:
  std::chrono::time_point<std::chrono::system_clock> _starting_point;
  std::chrono::time_point<std::chrono::system_clock> _end_point;
  bool _runing = false;
};

class Timer {
public:
  Timer() { SetStartingPointNow(); }
  void SetStartingPointNow() {
    _starting_point = clock();
    _runing = true;
  }
  void StopTimer() {
    if (_runing) {
      _end_point = clock();
      _runing = false;
    }
  }
  double ElapsedSeconds() {
    clock_t end_point;
    if (_runing) {
      end_point = clock();
    } else {
      end_point = _end_point;
    }
    return static_cast<double>(end_point - _starting_point) / CLOCKS_PER_SEC;
  }

private:
  clock_t _starting_point;
  clock_t _end_point;
  bool _runing = false;
};

//////////////////////////////////////////////////////////////////////////
// Implementation
//////////////////////////////////////////////////////////////////////////
template <typename T>
vector<T> ReadData(std::string filename, std::string input_type, unsigned n,
                   unsigned line_offset) {
  ifstream ifs(filename);
  vector<T> result;
  if (!ifs.is_open()) {
    std::cerr << "Cannot open file! (filename: " << filename << ")\n";
    exit(-1);
  }
  unsigned org_n = n;
  if (n == 0)
    n = std::numeric_limits<unsigned>::max();
  unsigned count = 0;
  T x = 0;
  if (input_type == "simple") {
    while (count < n + line_offset && ifs >> x) {
      if (count >= line_offset)
        result.push_back(x);
      ++count;
    }
  } else if (input_type == "multirecord") {
    std::string line;
    while (count < n + line_offset && std::getline(ifs, line)) {
      std::istringstream iss(line);
      if (!(iss >> x >> x)) {
        std::cerr << "Input file foramt error. \n";
        exit(-1);
      }
      if (count >= line_offset)
        result.push_back(x);
      ++count;
    }
  } else {
    MSG_ERROR(-1, "Invalid input-type: %s\n", input_type.c_str());
  }
  if (org_n && org_n != result.size()) {
    MSG_ERROR(-1, "Cannot read %u element from file %s. Only %zu read.\n",
              org_n, filename.c_str(), result.size());
  }
  ifs.close();
  return result;
}

template <typename T, unsigned K>
vector<KDPoint<T, K>> GetKDPoints(typename vector<T>::const_iterator first,
                                  typename vector<T>::const_iterator last,
                                  int count) {
  const size_t n = last - first;

  if (n <= K) {
    std::cerr << "GetKDPoints(): Data length is too short (n = " << n;
    std::cerr << ", K = " << K << ")" << std::endl;
    return vector<KDPoint<T, K>>();
  } else {
    vector<KDPoint<T, K>> points(n - K + 1);
    for (size_t i = 0; i < n - K + 1; i++) {
      points[i] = KDPoint<T, K>(first + i, first + i + K, count);
    }
    return points;
  }
}
template <typename T, unsigned K>
vector<vector<KDPoint<T, K>>>
GetKDPointsSample(typename vector<T>::const_iterator first,
                  typename vector<T>::const_iterator last,
                  const vector<vector<unsigned>> &indices, int count,
                  bool presort, OutputLevel output_level) {
  const size_t n = last - first;

  if (n < K) {
    std::cerr << "GetKDPoints(): Data length is too short (n = " << n;
    std::cerr << ", K = " << K << ")" << std::endl;
    exit(-1);
  } else {
    const unsigned sample_num = indices.size();
    vector<vector<KDPoint<T, K>>> points(sample_num);
    vector<KDPoint<T, K>> orig_points = GetKDPoints<T, K>(first, last, count);
    vector<unsigned> rank2index(orig_points.size());
    for (size_t i = 0; i < orig_points.size(); i++)
      rank2index.at(i) = i;

    if (presort) {
      Timer timer;
      std::sort(rank2index.begin(), rank2index.end(),
                [&orig_points](unsigned i1, unsigned i2) {
                  return (orig_points[i1] < orig_points[i2]);
                });
      timer.StopTimer();
      if (output_level == Debug) {
        MSG_INFO("Time consumed in presort: %f\n", timer.ElapsedSeconds());
      }
    }

    for (unsigned i_index = 0; i_index < sample_num; ++i_index) {
      const unsigned sample_size = indices[i_index].size();
      points[i_index].resize(sample_size);
      for (unsigned i = 0; i < sample_size; i++) {
        unsigned index = indices[i_index][i];
        assert(index < n - K + 1 &&
               "Sample index exceeds the number of points");
        points[i_index][i] = orig_points[rank2index[index]];
      }
    }
    return points;
  }
}

template <typename T> T ComputeSum(const vector<T> &data) {
  unsigned n0 = 1 << 4;
  unsigned p = data.size() / n0;
  T sum = std::accumulate(data.cbegin() + p * n0, data.cend(), (T)0);
  if (p == 0)
    return sum;
  else {
    vector<T> temp_sum(p, 0);
    for (unsigned i = 0; i < p; i++) {
      temp_sum[i] = std::accumulate(data.cbegin() + i * n0,
                                    data.cbegin() + (i + 1) * n0, (T)0);
    }
    temp_sum.push_back(sum);
    return ComputeSum(temp_sum);
  }
}

template <typename T> double ComputeVariance(const vector<T> &data) {
  vector<long double> data_(data.cbegin(), data.cend());
  long double avg = ComputeSum(data_) / data.size();
  std::for_each(data_.begin(), data_.end(), [avg](long double &x) {
    x -= avg;
    x *= x;
  });
  long double sum = ComputeSum(data_);
#ifdef DEBUG
  std::cout << "avg: " << avg << std::endl;
  std::cout << "variance * size: " << sum << std::endl;
#endif
  sum /= data.size();
  return (double)sum;
}

template <typename T, unsigned K>
void CloseAuxiliaryPoints(vector<KDPoint<T, K>> &points,
                          const vector<unsigned> &rank2index) {
  const unsigned n = points.size();

  for (unsigned i = 0; i < n; i++) {
    if (rank2index[i] >= n - K + 1)
      points[i].set_count(0);
    else
      points[i].set_count(1);
  }
}

template <typename T, unsigned K>
void MergeRepeatedPoints(vector<KDPoint<T, K>> &points,
                         const vector<unsigned> &rank2index) {
  const unsigned n = points.size();

  unsigned k = 0, i = 0;
  ;
  while (i < n) {
    int count_auxiliary = 0;
    int count = 0;
    while (k < n && points[k] == points[i]) {
      // Auxiliary points.
      if (rank2index[k] >= n - K + 1)
        count_auxiliary++;
      count += points[k].count();
      points[k].set_count(0);
      k++;
    }
    points[k - 1].set_count(count - count_auxiliary);
    i = k;
  }
}

template <typename T, unsigned K>
vector<KDPoint<unsigned, K - 1>> Map2Grid(const vector<KDPoint<T, K>> &points,
                                          const vector<unsigned> &rank2index,
                                          bool skip_nocount) {
  const unsigned n = points.size();

  assert(n == rank2index.size());
  vector<KDPoint<unsigned, K - 1>> result(n);

  vector<unsigned> index2rank = GetInverseMap(rank2index);
  // Mapping q in the paper.
  vector<unsigned> rank2next(n + 1);
  for (unsigned i = 0; i < n; i++) {
    rank2next[i] = index2rank[(rank2index[i] + 1) % n];
  }
  if (skip_nocount) {
    for (unsigned i = 0; i < n; i++) {
      if (points[i].count() == 0)
        continue;
      result[i].set_count(points[i].count());
      unsigned grid = i;
      for (unsigned j = 0; j < K - 1; j++) {
        grid = rank2next[grid];
        result[i][j] = grid;
      }
    }
  } else {
    for (unsigned i = 0; i < n; i++) {
      result[i].set_count(points[i].count());
      unsigned grid = i;
      for (unsigned j = 0; j < K - 1; j++) {
        grid = rank2next[grid];
        result[i][j] = grid;
      }
    }
  }

  return result;
}

template <typename T, unsigned K>
Bounds GetRankBounds(const vector<KDPoint<T, K>> &points, T r) {
  size_t n = points.size();
  Bounds bounds(n);
  vector<T> data(n);
  for (size_t i = 0; i < n; i++)
    data[i] = points[i][0];

  size_t k = 0;
  for (size_t i = 0; i < n; i++) {
    while (data[k] + r < data[i])
      k++;
    bounds.lower_bounds[i] = k;
  }
  k = n - 1;
  for (size_t i = n; i > 0; i--) {
    while (data[k] - r > data[i - 1])
      k--;
    bounds.upper_bounds[i - 1] = k;
  }
  return bounds;
}

template <unsigned K>
Range<unsigned, K> GetHyperCube(const KDPoint<unsigned, K> &point,
                                const Bounds &bounds) {
  Range<unsigned, K> result;
  for (size_t i = 0; i < K; ++i) {
    result.lower_ranges[i] = bounds.lower_bounds[point[i]];
    result.upper_ranges[i] = bounds.upper_bounds[point[i]];
  }
  return result;
}

} // namespace sampen

#endif // !__UTILS__
