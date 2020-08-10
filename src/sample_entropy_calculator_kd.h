#ifndef __SAMPLE_ENTROPY_CALCULATOR_KD__
#define __SAMPLE_ENTROPY_CALCULATOR_KD__
#include <iostream>
#include <vector>

#include <algorithm>

#include "kdtree.h"
#include "utils.h"


namespace kdtree_mddc
{
using std::vector;

// The kd tree of Mao Dong's version. 
template<typename T, unsigned K>
class MatchedPairsCalculatorMao
{
public:
    MatchedPairsCalculatorMao(OutputLevel output_level) 
        : _output_level(output_level) {}
    long long ComputeA(typename vector<T>::const_iterator first,
                       typename vector<T>::const_iterator last,
                       T r);
private:
    OutputLevel _output_level;
};


template<typename T, unsigned K>
class SampleEntropyCalculatorMao : SampleEntropyCalculator<T, K>
{
public:
    SampleEntropyCalculatorMao(OutputLevel output_level) 
        : SampleEntropyCalculator<T, K>(output_level) {}
    double ComputeSampleEntropy(typename vector<T>::const_iterator first,
                                typename vector<T>::const_iterator last, T r)
    {
        const unsigned n = last - first;
        if (n <= K)
        {
            std::cerr << "Data length is too short (n = " << n;
            std::cerr << ", K = " << K << ")" << std::endl;
            return -1;
        }

        MatchedPairsCalculatorMao<T, K> b_cal(this->_output_level);
        MatchedPairsCalculatorMao<T, K + 1> a_cal(this->_output_level);
        double b = static_cast<double>(b_cal.ComputeA(first, last - 1, r));
        double a = static_cast<double>(a_cal.ComputeA(first, last, r));
        double result = ComputeSampen(a, b, n - K, K, this->_output_level);
        return result;
    }
};

template<typename T, unsigned K>
class ABCalculatorLiu
{
public:
    ABCalculatorLiu(OutputLevel output_level) : _output_level(output_level) {}
    vector<long long> ComputeAB(typename vector<T>::const_iterator first, 
                                typename vector<T>::const_iterator last,
                                T r);
private:
    OutputLevel _output_level;
};

/*
* @brief This class calculates matched pairs like the class MatchedPairCalculatorMao,
* except that only one kd tree is used.
*/
template<typename T, unsigned K>
class SampleEntropyCalculatorLiu : public SampleEntropyCalculator<T, K>
{
public:
    SampleEntropyCalculatorLiu(OutputLevel output_level) 
        : SampleEntropyCalculator<T, K>(output_level) {}
    double ComputeSampleEntropy(typename vector<T>::const_iterator first,
                                typename vector<T>::const_iterator last, T r)
    {
        const unsigned n = last - first;
        if (n <= K)
        {
            std::cerr << "Data length is too short (n = " << n;
            std::cerr << ", K = " << K << ")" << std::endl;
            return -1;
        }
        ABCalculatorLiu<T, K> abc(this->_output_level);
        vector<long long> result = abc.ComputeAB(first, last, r);
        double entropy = ComputeSampen(result[0], result[1], 
                                       n - K, K, this->_output_level);
        return entropy;
    }
};


//////////////////////////////////////////////////////////////////////////
// Implementations
//////////////////////////////////////////////////////////////////////////

template<typename T, unsigned K>
long long MatchedPairsCalculatorMao<T, K>::ComputeA(
    typename vector<T>::const_iterator first,
    typename vector<T>::const_iterator last, T r)
{
    const size_t n = last - first;
    vector<T> data_(first, last);
    // Add K - 1 auxiliary points. 
    T minimum = *std::min_element(first, last);
    // TODO: !! For debug here. 
    for (size_t i = 0; i < K - 1; i++) data_.push_back(minimum);
    // Construct Points and merge repeated points. 
    const vector<KDPoint<T, K> > points = GetKDPoints<T, K>(data_);

    vector<KDPoint<T, K> > sorted_points(points);
    // The mapping p, from rank to original index 
    vector<unsigned> rank2index(n);
    for (size_t i = 0; i < n; i++) rank2index.at(i) = i;
    std::sort(rank2index.begin(), rank2index.end(),
              [&points](unsigned i1, unsigned i2) { return (points[i1] < points[i2]); });
    for (size_t i = 0; i < n; i++) sorted_points[i] = points[rank2index[i]];

    MergeRepeatedPoints(sorted_points, rank2index);

    const Bounds bounds = GetRankBounds(sorted_points, r);
    const vector<KDPoint<unsigned, K - 1> > points_grid = Map2Grid(sorted_points, rank2index);

    // Construct kd tree.
    vector<KDPoint<unsigned, K - 1> > points_count;
    vector<unsigned> points_count_indices;
    for (unsigned i = 0; i < n; i++)
    {
        if (points_grid[i].count())
        {
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

    for (unsigned i = 0; i < n_count - 1; i++)
    {
        const unsigned rank1 = points_count_indices[i];
        unsigned upperbound = bounds.upper_bounds[rank1];
        long long count_repeated = static_cast<long long>(points_count[i].count());
        result += (count_repeated - 1) * count_repeated / 2;

        if (upperbound < points_count_indices[i + 1]) continue;
        // Close current node. 
        tree.Close(i);

        // Update tree. 
        if (upperbound_prev < rank1) upperbound_prev = rank1;
        unsigned j = i + 1;
        while (j < n_count && points_count_indices[j] <= upperbound_prev) ++j;
        while (j < n_count && points_count_indices[j] <= upperbound)
        {
            tree.UpdateCount(j, points_count[j].count());
            ++num_opened;
            ++j;
        }

        const Range<unsigned, K - 1> range = GetHyperCube(points_count[i], bounds);
        result += tree.CountRange(range, num_nodes) * count_repeated;
        ++num_countrange_called;
        upperbound_prev = upperbound;
    }

    if (_output_level)
    {
        std::cout << "The number of nodes (K = " << K << "): ";
        std::cout << tree.num_nodes() << std::endl;;
        std::cout << "The number of leaf nodes (K = " << K << "): ";
        std::cout << n_count << std::endl;
        std::cout << "The number of calls for CountRange(): ";
        std::cout << num_countrange_called << std::endl;
        std::cout << "The number of times to open node: ";
        std::cout << num_opened << std::endl;
        std::cout << "The number of nodes visited (K = " << K << "): ";
        std::cout << num_nodes << std::endl;
    }
    return result;
}


template<typename T, unsigned K>
inline vector<long long> 
ABCalculatorLiu<T, K>::ComputeAB(typename vector<T>::const_iterator first, 
                                 typename vector<T>::const_iterator last, T r)
{
    const unsigned n = last - first;
    vector<T> data_(first, last);
    // Add K - 1 auxiliary points. 
    T minimum = *std::min_element(first, last);
    for (size_t i = 0; i < K; i++) data_.push_back(minimum);
    // Construct Points and merge repeated points. 
    const vector<KDPoint<T, K + 1> > points = GetKDPoints<T, K + 1>(data_);

    vector<KDPoint<T, K + 1> > sorted_points(points);
    // The mapping p, from rank to original index 
    vector<unsigned> rank2index(n);
    for (size_t i = 0; i < n; i++) rank2index.at(i) = i;
    std::sort(rank2index.begin(), rank2index.end(),
              [&points](unsigned i1, unsigned i2) { return (points[i1] < points[i2]); });
    for (size_t i = 0; i < n; i++) sorted_points[i] = points[rank2index[i]];

    //MergeRepeatedPoints(sorted_points, rank2index);
    CloseAuxiliaryPoints(sorted_points, rank2index);

    const Bounds bounds = GetRankBounds(sorted_points, r);
    const vector<KDPoint<unsigned, K> > points_grid = Map2Grid(sorted_points, rank2index);

    // Construct kd tree.
    vector<KDPoint<unsigned, K> > points_count;
    vector<unsigned> points_count_indices;
    for (unsigned i = 0; i < n; i++)
    {
        if (points_grid[i].count())
        {
            points_count.push_back(points_grid[i]);
            points_count_indices.push_back(i);
        }
    }
    KDTree2K<unsigned, K - 1> tree(points_count, _output_level);

    // Perform counting. 
    vector<long long> result({ 0, 0 });
    // The number of nodes has been visited. 
    long long num_nodes = 0;
    long long num_countrange_called = 0;
    long long num_opened = 0;
    unsigned upperbound_prev = 0;

    const unsigned n_count = points_count.size();

    for (unsigned i = 0; i < n_count - 1; i++)
    {
        const unsigned rank1 = points_count_indices[i];
        unsigned upperbound = bounds.upper_bounds[rank1];
        if (upperbound < points_count_indices[i + 1]) continue;
        // Close current node. 
        tree.Close(i);

        // Update tree. 
        if (upperbound_prev < rank1) upperbound_prev = rank1;
        unsigned j = i + 1;
        while (j < n_count && points_count_indices[j] <= upperbound_prev) ++j;
        while (j < n_count && points_count_indices[j] <= upperbound)
        {
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

    if (_output_level)
    {
        std::cout << "The number of nodes (K = " << K << "): ";
        std::cout << tree.num_nodes() << std::endl;;
        std::cout << "The number of leaf nodes (K = " << K << "): ";
        std::cout << n_count << std::endl;
        std::cout << "The number of calls for CountRange(): ";
        std::cout << num_countrange_called << std::endl;
        std::cout << "The number times to open node: ";
        std::cout << num_opened << std::endl;
        std::cout << "The number of nodes visited (K = " << K << "): ";
        std::cout << num_nodes << std::endl;
    }

    return result;
}

} // namespace kdtree_mddc

#endif // !__SAMPLE_ENTROPY_CALCULATOR_KD__
