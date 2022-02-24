//
// Created by Phree on 2021/12/9.
//

#include "sample_entropy_calculator_kd.h"

namespace sampen {

template <typename T, unsigned K>
long long MatchedPairsCalculatorMao<T, K>::ComputeA(
    typename vector<T>::const_iterator first,
    typename vector<T>::const_iterator last, T r) {
  const size_t n = last - first;
  vector<T> data_(first, last);
  // Add K - 1 auxiliary points.
  T minimum = *std::min_element(first, last);
  for (size_t i = 0; i < K - 1; i++)
    data_.push_back(minimum);
  // Construct Points and merge repeated points.
  vector<KDPoint<T, K>> points =
      GetKDPoints<T, K>(data_.cbegin(), data_.cend(), 1);
  for (size_t i = points.size() - K + 1; i < points.size(); ++i) {
    points[i].set_count(0);
  }
  vector<KDPoint<T, K>> sorted_points(points);
  // The mapping p, from rank to original index
  vector<unsigned> rank2index(n);
  for (size_t i = 0; i < n; i++)
    rank2index.at(i) = i;

  Timer timer;
  std::sort(rank2index.begin(), rank2index.end(),
            [&points](unsigned i1, unsigned i2) {
              return (points[i1] < points[i2]);
            });
  timer.StopTimer();
  if (_output_level >= Info) {
    std::cout << "[INFO] Time consumed in presorting: "
              << timer.ElapsedSeconds() << "s\n";
  }

  for (size_t i = 0; i < n; i++)
    sorted_points[i] = points[rank2index[i]];

  // MergeRepeatedPoints(sorted_points, rank2index);

  const Bounds bounds = GetRankBounds(sorted_points, r);
  const vector<KDPoint<unsigned, K - 1>> points_grid =
      Map2Grid(sorted_points, rank2index);

  // Construct kd tree.
  vector<KDPoint<unsigned, K - 1>> points_count;
  vector<unsigned> points_count_indices;
  for (unsigned i = 0; i < n; i++) {
    if (points_grid[i].count()) {
      points_count.push_back(points_grid[i]);
      points_count_indices.push_back(i);
    }
  }
  KDCountingTree2K<unsigned, K - 1> tree(points_count, _output_level);

  // Perform counting.
  long long result = 0;
  // The number of nodes has been visited.
  long long num_nodes = 0;
  long long num_countrange_called = 0;
  long long num_opened = 0;
  unsigned upperbound_prev = 0;

  const unsigned n_count = points_count.size();

  timer.SetStartingPointNow();
  for (unsigned i = 0; i < n_count - 1; i++) {
    // Close current node.
    tree.Close(i);

    const unsigned rank1 = points_count_indices[i];

    unsigned upperbound = bounds.upper_bounds[rank1];
    long long count_repeated = static_cast<long long>(points_count[i].count());
    result += (count_repeated - 1) * count_repeated / 2;

    if (upperbound < points_count_indices[i + 1])
      continue;

    // Update tree.
    if (upperbound_prev < rank1)
      upperbound_prev = rank1;
    unsigned j = i + 1;
    while (j < n_count && points_count_indices[j] <= upperbound_prev)
      ++j;
    while (j < n_count && points_count_indices[j] <= upperbound) {
      tree.UpdateCount(j, points_count[j].count());
      ++num_opened;
      ++j;
    }

    const Range<unsigned, K - 1> range = GetHyperCube(points_count[i], bounds);
    long long current_count =
        tree.CountRange(range, num_nodes) * count_repeated;
    result += current_count;
    ++num_countrange_called;
    upperbound_prev = upperbound;
  }
  timer.StopTimer();

  if (_output_level >= Info) {
    std::cout << "[INFO] Time consumed in range counting: "
              << timer.ElapsedSeconds() << " seconds\n";
  }
  if (_output_level == Debug) {
    std::cout << "[INFO] The number of nodes (K = " << K << "): ";
    std::cout << tree.num_nodes() << std::endl;
    std::cout << "[INFO] The number of leaf nodes (K = " << K << "): ";
    std::cout << n_count << std::endl;
    std::cout << "[INFO] The number of calls for CountRange(): ";
    std::cout << num_countrange_called << std::endl;
    std::cout << "[INFO] The number of times to open node: ";
    std::cout << num_opened << std::endl;
    std::cout << "[INFO] The number of nodes visited (K = " << K << "): ";
    std::cout << num_nodes << std::endl;
  }
  return result;
}

template <typename T, unsigned K>
long long MatchedPairsCalculatorSampling2<T, K>::ComputeA(
    typename vector<T>::const_iterator first,
    typename vector<T>::const_iterator last, T r,
    std::vector<unsigned> &sample_indices) {
  const size_t n = last - first;
  assert(n - K + 1 >= sample_indices.size());
  vector<T> data_(first, last);
  // Add K - 1 auxiliary points.
  T minimum = *std::min_element(first, last);
  for (size_t i = 0; i < K - 1; ++i)
    data_.push_back(minimum);

  // Construct points in k-dimensional space and merge repeated points.
  vector<KDPoint<T, K> > points =
      GetKDPoints<T, K>(data_.cbegin(), data_.cend(), 1);
  for (size_t i = points.size() - K + 1; i < points.size(); ++i) {
    points[i].set_count(0);
  }

  vector<KDPoint<T, K> > sorted_points(points);
  // The mapping p, from rank to original index
  vector<unsigned> rank2index(n);
  for (size_t i = 0; i < n; i++)
    rank2index.at(i) = i;
  Timer timer;
  std::sort(rank2index.begin(), rank2index.end(),
            [&points](unsigned i1, unsigned i2) {
              return (points[i1] < points[i2]);
            });
  timer.StopTimer();
  if (_output_level >= Info) {
    std::cout << "[INFO] Time consumed in presorting: "
              << timer.ElapsedSeconds() << "s\n";
  }

  for (size_t i = 0; i < n; i++)
    sorted_points[i] = points[rank2index[i]];

  const Bounds bounds = GetRankBounds(sorted_points, r);
  // Map values of each coordinate to the rank given by sorting.
  // Since the value at first dimension equal to the index of that point
  // in the sorted array, we can reduce the dimension of the points by 1.
  const vector<KDPoint<unsigned, K - 1>> points_grid =
      Map2Grid(sorted_points, rank2index);

  // Construct kd tree.
  vector<KDPoint<unsigned, K - 1>> points_count;
  vector<unsigned> points_count_indices;
  // Only use points with positive count to construct the kd tree.
  for (unsigned i = 0; i < n; i++) {
    if (points_grid[i].count()) {
      points_count.push_back(points_grid[i]);
      points_count_indices.push_back(i);
    }
  }
  KDCountingTree2K<unsigned, K - 1> tree(points_count, _output_level);

  // Perform counting.
  long long result = 0;
  // The number of nodes has been visited.
  long long num_nodes = 0;
  long long num_countrange_called = 0;
  long long num_opened = 0;
  unsigned upperbound_prev = 0;

  const unsigned n_count = points_count.size();
  assert(n_count == n - K + 1);
  int i_sample_indices = 0;
  std::sort(sample_indices.begin(), sample_indices.end());
  for (size_t i = 0; i < sample_indices.size() - 1; ++i) {
    assert(sample_indices[i + 1] > sample_indices[i]);
  }
  timer.SetStartingPointNow();
  for (unsigned i = 0; i < n_count; ++i) {
    const unsigned rank1 = points_count_indices[i];
    unsigned upperbound = bounds.upper_bounds[rank1];

    // Close nodes whose value of the first dimension are outside bounds.
    if (i_sample_indices > 0) {
      for (unsigned k = sample_indices[i_sample_indices - 1] + 1; k <= i; ++k) {
        tree.Close(k);
      }
    }

    if (sample_indices[i_sample_indices] != i) {
      continue;
    }
    if (upperbound < points_count_indices[i + 1]) {
      ++i_sample_indices;
      continue;
    }

    // Update tree.
    if (upperbound_prev < rank1)
      upperbound_prev = rank1;
    unsigned j = i + 1;
    while (j < n_count && points_count_indices[j] <= upperbound_prev)
      ++j;
    while (j < n_count && points_count_indices[j] <= upperbound) {
      tree.UpdateCount(j, points_count[j].count());
      ++num_opened;
      ++j;
    }

    const Range<unsigned, K - 1> range = GetHyperCube(points_count[i], bounds);
    result += tree.CountRange(range, num_nodes);
    ++num_countrange_called;
    upperbound_prev = upperbound;
    ++i_sample_indices;
  }
  timer.StopTimer();

  assert(i_sample_indices == static_cast<int>(sample_indices.size()));

  if (_output_level >= Info) {
    std::cout << "[INFO] Time consumed in range counting: "
              << timer.ElapsedSeconds() << " seconds\n";
  }

  if (_output_level == Debug) {
    std::cout << "[DEBUG] The number of nodes (K = " << K << "): ";
    std::cout << tree.num_nodes() << std::endl;
    std::cout << "[DEBUG] The number of leaf nodes (K = " << K << "): ";
    std::cout << n_count << std::endl;
    std::cout << "[DEBUG] The number of calls for CountRange(): ";
    std::cout << num_countrange_called << std::endl;
    std::cout << "[DEBUG] The number of times to open node: ";
    std::cout << num_opened << std::endl;
    std::cout << "[DEBUG] The number of nodes visited (K = " << K << "): ";
    std::cout << num_nodes << std::endl;
  }
  return result;
}
vector<unsigned>
MergeRepeatedIndices(typename vector<unsigned>::const_iterator first,
                     typename vector<unsigned>::const_iterator last) {
  size_t length = last - first;
  vector<unsigned> counts(length, 0);
  unsigned i = 0, k = 0;
  while (i < length) {
    unsigned count = 0;
    while (k < length && *(first + k) == *(first + i)) {
      ++k;
      ++count;
    }
    counts[i] = count;
    i = k;
  }
  return counts;
}

template <typename T, unsigned K>
vector<long long> MatchedPairsCalculatorSampling<T, K>::ComputeA(
    typename vector<T>::const_iterator first,
    typename vector<T>::const_iterator last, unsigned sample_num,
    const vector<unsigned> &indices, T r) {
  const size_t n = last - first;
  const size_t sample_size = indices.size() / sample_num;
  vector<T> data_(first, last);
  // Add K - 1 auxiliary points.
  T minimum = *std::min_element(first, last);
  // TODO: !! For debug here.
  for (size_t i = 0; i < K - 1; i++)
    data_.push_back(minimum);
  // Construct Points and merge repeated points.
  vector<KDPoint<T, K>> points =
      GetKDPoints<T, K>(data_.cbegin(), data_.cend(), 0);
  // The mapping p, from rank to original index
  vector<unsigned> rank2index(n);
  for (size_t i = 0; i < n; i++)
    rank2index.at(i) = i;

  Timer timer;
  std::sort(rank2index.begin(), rank2index.end(),
            [&points](unsigned i1, unsigned i2) {
              return (points[i1] < points[i2]);
            });
  timer.StopTimer();
  if (_output_level == Info) {
    std::cout << "[INFO] Time consumed in presorting: "
              << timer.ElapsedSeconds() << "s\n";
  }

  vector<KDPoint<T, K>> sorted_points(points.cbegin(), points.cend());
  for (size_t i = 0; i < n; i++)
    sorted_points[i] = points[rank2index[i]];

  // MergeRepeatedPoints(sorted_points, rank2index);

  const Bounds bounds = GetRankBounds(sorted_points, r);
  vector<KDPoint<unsigned, K - 1>> points_grid =
      Map2Grid(sorted_points, rank2index, false);

  vector<vector<unsigned>> sample_groups(n);
  for (unsigned i = 0; i < sample_num * sample_size; ++i) {
    unsigned index = indices[i];
    points_grid[index].increase_count(1);
    sample_groups[index].push_back(i / sample_size);
  }

  // Construct kd tree.
  vector<KDPoint<unsigned, K - 1>> points_count;
  vector<unsigned> points_count_indices;
  vector<vector<unsigned>> indices_count(sample_num);
  for (unsigned i = 0; i < n; i++) {
    if (points_grid[i].count()) {
      for (unsigned j = 0; j < sample_groups[i].size(); ++j) {
        indices_count[sample_groups[i][j]].push_back(points_count.size());
      }
      points_count.push_back(points_grid[i]);
      points_count_indices.push_back(i);
    }
  }

  KDCountingTree2K<unsigned, K - 1> tree(points_count, _output_level);

  // Perform counting.
  // The number of nodes has been visited.
  long long num_nodes = 0;
  long long num_countrange_called = 0;
  long long num_opened = 0;
  unsigned upperbound_prev = 0;

  const unsigned n_count = points_count.size();

  timer.SetStartingPointNow();
  for (unsigned i = 0; i < n_count; ++i) {
    points_count[i].set_count(0);
  }
  vector<long long> results(sample_num);
  for (unsigned k = 0; k < sample_num; ++k) {
    for (unsigned i = 0; i < indices_count[k].size(); ++i) {
      points_count[indices_count[k][i]].increase_count(1);
    }
    long long result = 0;
    for (unsigned i = 0; i < n_count - 1; ++i) {
      if (points_count[i].count() == 0)
        continue;

      // Close current node.
      tree.Close(i);

      const unsigned rank1 = points_count_indices[i];
      unsigned upperbound = bounds.upper_bounds[rank1];
      long long count_repeated =
          static_cast<long long>(points_count[i].count());
      result += (count_repeated - 1) * count_repeated / 2;

      if (upperbound < points_count_indices[i + 1])
        continue;

      // Update tree.
      if (upperbound_prev < rank1)
        upperbound_prev = rank1;
      unsigned j = i + 1;
      while (j < n_count && points_count_indices[j] <= upperbound_prev)
        ++j;
      while (j < n_count && points_count_indices[j] <= upperbound) {
        if (points_count[j].count()) {
          tree.UpdateCount(j, points_count[j].count());
          ++num_opened;
        }
        ++j;
      }

      const Range<unsigned, K - 1> range =
          GetHyperCube(points_count[i], bounds);
      result += tree.CountRange(range, num_nodes) * count_repeated;
      ++num_countrange_called;
      upperbound_prev = upperbound;
    }
    tree.Close(n_count - 1);
    results[k] = result;

    // Close points for current sample group.
    for (unsigned i = 0; i < indices_count[k].size(); ++i) {
      points_count[indices_count[k][i]].set_count(0);
    }
  }
  timer.StopTimer();

  if (_output_level >= Info) {
    std::cout << "[INFO] Time consumed in range counting: "
              << timer.ElapsedSeconds() << " seconds\n";
  }
  if (_output_level == Debug) {
    std::cout << "[DEBUG] The number of nodes (K = " << K << "): ";
    std::cout << tree.num_nodes() << std::endl;
    std::cout << "[DEBUG] The number of leaf nodes (K = " << K << "): ";
    std::cout << n_count << std::endl;
    std::cout << "[DEBUG] The number of calls for CountRange(): ";
    std::cout << num_countrange_called << std::endl;
    std::cout << "[DEBUG] The number of times to open node: ";
    std::cout << num_opened << std::endl;
    std::cout << "[DEBUG] The number of nodes visited (K = " << K << "): ";
    std::cout << num_nodes << std::endl;
    std::cout << "[DEBUG] The numbers of the matched pairs: \n";
    for (unsigned i = 0; i < sample_num; ++i) {
      if (i)
        std::cout << ", " << results[i];
      else
        std::cout << results[i];
    }
    std::cout << std::endl;
    for (unsigned i = 0; i < sample_num; ++i) {
      std::cout << "[DEBUG] Sample indices (count): \n";
      for (unsigned j = 0; j < sample_size; ++j) {
        std::cout << indices_count[i][j] << ", ";
      }
      std::cout << std::endl;
    }
  }
  return results;
}

template <typename T, unsigned K>
inline vector<long long>
ABCalculatorLiu<T, K>::ComputeAB(typename vector<T>::const_iterator first,
                                 typename vector<T>::const_iterator last, T r) {
  const unsigned n = last - first;
  vector<T> data_(first, last);
  // Add K - 1 auxiliary points.
  T minimum = *std::min_element(first, last);
  for (size_t i = 0; i < K; i++)
    data_.push_back(minimum);
  // Construct Points and merge repeated points.
  const vector<KDPoint<T, K + 1>> points =
      GetKDPoints<T, K + 1>(data_.cbegin(), data_.cend());

  vector<KDPoint<T, K + 1>> sorted_points(points);
  // The mapping p, from rank to original index
  vector<unsigned> rank2index(n);
  for (size_t i = 0; i < n; i++)
    rank2index.at(i) = i;

  Timer timer;
  std::sort(rank2index.begin(), rank2index.end(),
            [&points](unsigned i1, unsigned i2) {
              return (points[i1] < points[i2]);
            });
  for (size_t i = 0; i < n; i++)
    sorted_points[i] = points[rank2index[i]];
  timer.StopTimer();
  if (_output_level >= Info) {
    std::cout << "[INFO] Time consumed in presorting: "
              << timer.ElapsedSeconds() << " seconds\n";
  }

  // MergeRepeatedPoints(sorted_points, rank2index);
  CloseAuxiliaryPoints(sorted_points, rank2index);

  const Bounds bounds = GetRankBounds(sorted_points, r);
  const vector<KDPoint<unsigned, K>> points_grid =
      Map2Grid(sorted_points, rank2index);

  // Construct kd tree.
  vector<KDPoint<unsigned, K>> points_count;
  vector<unsigned> points_count_indices;
  for (unsigned i = 0; i < n; i++) {
    if (points_grid[i].count()) {
      points_count.push_back(points_grid[i]);
      points_count_indices.push_back(i);
    }
  }
  KDTree2K<unsigned, K - 1> tree(points_count, _output_level);

  // Perform counting.
  vector<long long> result({0, 0});
  // The number of nodes has been visited.
  long long num_nodes = 0;
  long long num_countrange_called = 0;
  long long num_opened = 0;
  unsigned upperbound_prev = 0;

  const unsigned n_count = points_count.size();

  timer.SetStartingPointNow();
  for (unsigned i = 0; i < n_count - 1; i++) {
    // Close current node.
    tree.Close(i);

    const unsigned rank1 = points_count_indices[i];
    unsigned upperbound = bounds.upper_bounds[rank1];

    if (upperbound < points_count_indices[i + 1])
      continue;
    // Update tree.
    if (upperbound_prev < rank1)
      upperbound_prev = rank1;
    unsigned j = i + 1;
    while (j < n_count && points_count_indices[j] <= upperbound_prev)
      ++j;
    while (j < n_count && points_count_indices[j] <= upperbound) {
      tree.UpdateCount(j, points_count[j].count());
      ++num_opened;
      ++j;
    }

    const Range<unsigned, K> range = GetHyperCube(points_count[i], bounds);
    vector<long long> ab = tree.CountRange(range, num_nodes);

    result[0] += ab[0];
    result[1] += ab[1];
    ++num_countrange_called;
    upperbound_prev = upperbound;
  }
  timer.StopTimer();

  if (_output_level >= Info) {
    std::cout << "[INFO] Time consumed in range counting: "
              << timer.ElapsedSeconds() << " seconds\n";
  }
  if (_output_level == Debug) {
    std::cout << "[DEBUG] The number of nodes (K = " << K << "): ";
    std::cout << tree.num_nodes() << std::endl;
    std::cout << "[DEBUG] The number of leaf nodes (K = " << K << "): ";
    std::cout << n_count << std::endl;
    std::cout << "[DEBUG] The number of calls for CountRange(): ";
    std::cout << num_countrange_called << std::endl;
    std::cout << "[DEBUG] The number times to open node: ";
    std::cout << num_opened << std::endl;
    std::cout << "[DEBUG] The number of nodes visited (K = " << K << "): ";
    std::cout << num_nodes << std::endl;
  }

  return result;
}


template <typename T, unsigned K>
inline vector<long long>
ABCalculatorRKD<T, K>::ComputeAB(typename vector<T>::const_iterator first,
                                 typename vector<T>::const_iterator last, T r) {
  const unsigned n = last - first;
  vector<T> data_(first, last);
  // Add K - 1 auxiliary points.
  T minimum = *std::min_element(first, last);
  for (size_t i = 0; i < K; i++)
    data_.push_back(minimum);
  // Construct Points and merge repeated points.
  const vector<KDPoint<T, K + 1>> points =
      GetKDPoints<T, K + 1>(data_.cbegin(), data_.cend());

  vector<KDPoint<T, K + 1>> sorted_points(points);
  // The mapping p, from rank to original index
  vector<unsigned> rank2index(n);
  for (size_t i = 0; i < n; i++)
    rank2index.at(i) = i;

  Timer timer;
  std::sort(rank2index.begin(), rank2index.end(),
            [&points](unsigned i1, unsigned i2) {
              return (points[i1] < points[i2]);
            });
  for (size_t i = 0; i < n; i++)
    sorted_points[i] = points[rank2index[i]];
  timer.StopTimer();
  if (_output_level >= Info) {
    std::cout << "[INFO] Time consumed in presorting: "
              << timer.ElapsedSeconds() << " seconds\n";
  }

  // MergeRepeatedPoints(sorted_points, rank2index);
  CloseAuxiliaryPoints(sorted_points, rank2index);

  const Bounds bounds = GetRankBounds(sorted_points, r);
  const vector<KDPoint<unsigned, K>> points_grid =
      Map2Grid(sorted_points, rank2index);

  // Construct kd tree.
  vector<KDPointRKD<unsigned, K> > points_count;
  vector<unsigned> points_count_indices;
  for (unsigned i = 0; i < n; i++) {
    if (points_grid[i].count()) {
      points_count.push_back(KDPointRKD<unsigned, K>(points_grid[i]));
      points_count_indices.push_back(i);
    }
  }
  RangeKDTree2K<unsigned, K - 1> tree(points_count, 0, _output_level);

  // Perform counting.
  vector<long long> result({0, 0});
  // The number of nodes has been visited.
  long long num_nodes = 0;
  long long num_countrange_called = 0;
  long long num_opened = 0;
  unsigned upperbound_prev = 0;

  const unsigned n_count = points_count.size();

  timer.SetStartingPointNow();
  for (unsigned i = 0; i < n_count - 1; i++) {
    // Close current node.
    tree.Close(i);

    const unsigned rank1 = points_count_indices[i];
    unsigned upperbound = bounds.upper_bounds[rank1];

    if (upperbound < points_count_indices[i + 1])
      continue;
    // Update tree.
    if (upperbound_prev < rank1)
      upperbound_prev = rank1;
    unsigned j = i + 1;
    while (j < n_count && points_count_indices[j] <= upperbound_prev)
      ++j;
    while (j < n_count && points_count_indices[j] <= upperbound) {
      tree.UpdateCount(j, points_count[j].count());
      ++num_opened;
      ++j;
    }

    const Range<unsigned, K> range = GetHyperCube(points_count[i], bounds);
    vector<long long> ab = tree.CountRange(range, num_nodes);

    result[0] += ab[0];
    result[1] += ab[1];
    ++num_countrange_called;
    upperbound_prev = upperbound;
  }
  timer.StopTimer();

  if (_output_level >= Info) {
    std::cout << "[INFO] Time consumed in range counting: "
              << timer.ElapsedSeconds() << " seconds\n";
  }
  if (_output_level == Debug) {
    std::cout << "[DEBUG] The number of nodes (K = " << K << "): ";
    std::cout << tree.num_nodes() << std::endl;
    std::cout << "[DEBUG] The number of leaf nodes (K = " << K << "): ";
    std::cout << n_count << std::endl;
    std::cout << "[DEBUG] The number of calls for CountRange(): ";
    std::cout << num_countrange_called << std::endl;
    std::cout << "[DEBUG] The number times to open node: ";
    std::cout << num_opened << std::endl;
    std::cout << "[DEBUG] The number of nodes visited (K = " << K << "): ";
    std::cout << num_nodes << std::endl;
  }

  return result;
}


template <typename T, unsigned K>
vector<long long> ABCalculatorSamplingLiu<T, K>::ComputeAB(
    typename vector<T>::const_iterator first,
    typename vector<T>::const_iterator last, unsigned sample_num,
    const vector<unsigned> &indices, T r) {
  const unsigned n = last - first;
  vector<T> data_(first, last);
  const unsigned sample_size = indices.size() / sample_num;
  // Add K auxiliary points.
  T minimum = *std::min_element(first, last);
  for (size_t i = 0; i < K; i++)
    data_.push_back(minimum);
  // Construct Points and merge repeated points.
  const vector<KDPoint<T, K + 1>> points =
      GetKDPoints<T, K + 1>(data_.cbegin(), data_.cend(), 0);

  vector<KDPoint<T, K + 1>> sorted_points(points);
  // The mapping p, from rank to original index
  vector<unsigned> rank2index(n);
  for (size_t i = 0; i < n; i++)
    rank2index.at(i) = i;

  Timer timer;
  std::sort(rank2index.begin(), rank2index.end(),
            [&points](unsigned i1, unsigned i2) {
              return (points[i1] < points[i2]);
            });
  for (size_t i = 0; i < n; i++)
    sorted_points[i] = points[rank2index[i]];
  timer.StopTimer();
  if (_output_level >= Info) {
    std::cout << "[INFO] Time consumed in presorting: "
              << timer.ElapsedSeconds() << "s\n";
  }

  const Bounds bounds = GetRankBounds(sorted_points, r);
  vector<KDPoint<unsigned, K>> points_grid =
      Map2Grid(sorted_points, rank2index, false);

  vector<vector<unsigned>> sample_groups(n);
  for (unsigned i = 0; i < sample_num * sample_size; ++i) {
    unsigned index = indices[i];
    points_grid[index].increase_count(1);
    sample_groups[index].push_back(i / sample_size);
  }

  // Construct kd tree.
  vector<KDPoint<unsigned, K>> points_count;
  vector<unsigned> points_count_indices;
  vector<vector<unsigned>> indices_count(sample_num);
  for (unsigned i = 0; i < n; i++) {
    if (points_grid[i].count()) {
      for (unsigned j = 0; j < sample_groups[i].size(); ++j) {
        indices_count[sample_groups[i][j]].push_back(points_count.size());
      }
      points_count.push_back(points_grid[i]);
      points_count_indices.push_back(i);
    }
  }
  KDTree2K<unsigned, K - 1> tree(points_count, _output_level);

  // Perform counting.
  vector<long long> results(2 * sample_num);
  // The number of nodes has been visited.
  long long num_nodes = 0;
  long long num_countrange_called = 0;
  long long num_opened = 0;
  unsigned upperbound_prev = 0;

  const unsigned n_count = points_count.size();

  timer.SetStartingPointNow();
  for (unsigned i = 0; i < n_count; ++i) {
    points_count[i].set_count(0);
  }
  for (unsigned k = 0; k < sample_num; ++k) {
    for (unsigned i = 0; i < indices_count[k].size(); ++i) {
      points_count[indices_count[k][i]].increase_count(1);
    }
    long long result_a = 0;
    long long result_b = 0;
    for (unsigned i = 0; i < n_count - 1; ++i) {
      if (points_count[i].count() == 0)
        continue;

      // Close current node.
      tree.Close(i);

      const unsigned rank1 = points_count_indices[i];
      unsigned upperbound = bounds.upper_bounds[rank1];
      long long count_repeated =
          static_cast<long long>(points_count[i].count());
      result_a += (count_repeated - 1) * count_repeated / 2;
      result_b += (count_repeated - 1) * count_repeated / 2;

      if (upperbound < points_count_indices[i + 1])
        continue;

      // Update tree.
      if (upperbound_prev < rank1)
        upperbound_prev = rank1;
      unsigned j = i + 1;
      while (j < n_count && points_count_indices[j] <= upperbound_prev)
        ++j;
      while (j < n_count && points_count_indices[j] <= upperbound) {
        if (points_count[j].count()) {
          tree.UpdateCount(j, points_count[j].count());
          ++num_opened;
        }
        ++j;
      }

      const Range<unsigned, K> range = GetHyperCube(points_count[i], bounds);
      vector<long long> ab = tree.CountRange(range, num_nodes);

      result_a += ab[0] * count_repeated;
      result_b += ab[1] * count_repeated;
      ++num_countrange_called;
      upperbound_prev = upperbound;
    }
    results[2 * k] = result_a;
    results[2 * k + 1] = result_b;
    tree.Close(n_count - 1);

    for (unsigned i = 0; i < indices_count[k].size(); ++i) {
      points_count[indices_count[k][i]].set_count(0);
    }
  }
  timer.StopTimer();

  if (_output_level >= Info) {
    std::cout << "[INFO] Time consumed in range counting: "
              << timer.ElapsedSeconds() << " seconds\n";
  }
  if (_output_level == Debug) {
    std::cout << "[DEBUG] The number of nodes (K = " << K << "): ";
    std::cout << tree.num_nodes() << std::endl;
    std::cout << "[DEBUG] The number of leaf nodes (K = " << K << "): ";
    std::cout << n_count << std::endl;
    std::cout << "[DEBUG] The number of calls for CountRange(): ";
    std::cout << num_countrange_called << std::endl;
    std::cout << "[DEBUG] The number times to open node: ";
    std::cout << num_opened << std::endl;
    std::cout << "[DEBUG] The number of nodes visited (K = " << K << "): ";
    std::cout << num_nodes << std::endl;
  }

  return results;
}

#define INSTANTIATE_SAMPLE_ENTROPY_CALCULATOR(K) \
  INSTANTIATE_SAMPLE_ENTROPY_CALCULATOR_TYPE_K(int, K) \
  INSTANTIATE_SAMPLE_ENTROPY_CALCULATOR_TYPE_K(double, K)
  
#define INSTANTIATE_SAMPLE_ENTROPY_CALCULATOR_TYPE_K(TYPE, K) \
template class SampleEntropyCalculatorLiu<TYPE , K>; \
template class SampleEntropyCalculatorRKD<TYPE , K>; \
template class SampleEntropyCalculatorMao<TYPE , K>; \
template class SampleEntropyCalculatorSamplingLiu<TYPE , K>; \
template class SampleEntropyCalculatorSamplingMao<TYPE , K>; \
template class SampleEntropyCalculatorSamplingKDTree<TYPE , K>; \
template class MatchedPairsCalculatorMao<TYPE, K>; \
template class MatchedPairsCalculatorSampling<TYPE, K>; \
template class MatchedPairsCalculatorSampling2<TYPE, K>; \
template class ABCalculatorLiu<TYPE, K>; \
template class ABCalculatorRKD<TYPE, K>; \
template class ABCalculatorSamplingLiu<TYPE, K>;


INSTANTIATE_SAMPLE_ENTROPY_CALCULATOR(2);
INSTANTIATE_SAMPLE_ENTROPY_CALCULATOR(3);
INSTANTIATE_SAMPLE_ENTROPY_CALCULATOR(4);
INSTANTIATE_SAMPLE_ENTROPY_CALCULATOR(5);
INSTANTIATE_SAMPLE_ENTROPY_CALCULATOR(6);
INSTANTIATE_SAMPLE_ENTROPY_CALCULATOR(7);
INSTANTIATE_SAMPLE_ENTROPY_CALCULATOR(8);
INSTANTIATE_SAMPLE_ENTROPY_CALCULATOR(9);
INSTANTIATE_SAMPLE_ENTROPY_CALCULATOR(10);
INSTANTIATE_SAMPLE_ENTROPY_CALCULATOR(11);
} // namespace sampen