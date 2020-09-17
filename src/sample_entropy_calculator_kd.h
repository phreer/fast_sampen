#ifndef __SAMPLE_ENTROPY_CALCULATOR_KD__
#define __SAMPLE_ENTROPY_CALCULATOR_KD__
#include <iostream>
#include <vector>

#include <algorithm>
#include "random_sampler.h" 
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
class MatchedPairsCalculatorSample 
{
public: 
    MatchedPairsCalculatorSample(OutputLevel output_level) 
        : _output_level(output_level) {} 
    long long ComputeA(typename vector<T>::const_iterator first, 
                       typename vector<T>::const_iterator last, 
                       const vector<unsigned> &indices, 
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

        // Display 
        double norm = (n - K + 1) * (n - K); 
        double A_norm = static_cast<double>(a) / norm; 
        double B_norm = static_cast<double>(b) / norm; 
        std::cout << "A (norm): " << A_norm << ", B (norm): " 
            << B_norm << std::endl; 

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

        // Display 
        double norm = (n - K + 1) * (n - K); 
        double A_norm = static_cast<double>(result[0]) / norm; 
        double B_norm = static_cast<double>(result[1]) / norm; 
        std::cout << "A (norm): " << A_norm << ", B (norm): " 
            << B_norm << std::endl; 

        double entropy = ComputeSampen(result[0], result[1], 
                                       n - K, K, this->_output_level);
        return entropy;
    }
};

template<typename T, unsigned K>
class SampleEntropyCalculatorS : public SampleEntropyCalculator<T, K> 
{
public:
    SampleEntropyCalculatorS(unsigned sample_size, unsigned sample_num, 
                             uniform_int_generator::random_type rtype, 
                             bool random_, OutputLevel output_level) 
        : _sample_size(sample_size), _sample_num(sample_num), _rtype(rtype), 
        _random(random_), SampleEntropyCalculator<T, K>(output_level) {} 
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
        uniform_int_generator uig(0, n - 2, _rtype, _random); 
        MatchedPairsCalculatorSample<T, K> b_cal(this->_output_level);
        MatchedPairsCalculatorSample<T, K + 1> a_cal(this->_output_level);
        vector<double> ABs(2 * _sample_num); 

        for (unsigned i = 0; i < _sample_num; ++i) 
        {
            vector<unsigned> indices(_sample_size); 
            for (unsigned j = 0; j < _sample_size; ++j) 
            {
                indices[j] = uig.get(); 
            }
            ABs[2 * i] = a_cal.ComputeA(first, last - 1, indices, r); 
            ABs[2 * i + 1] = b_cal.ComputeA(first, last, indices, r); 
        }
        double A = 0, B = 0; 
        for (unsigned i = 0; i < _sample_num; ++i) 
        {
            A += ABs[2 * i]; 
            B += ABs[2 * i + 1]; 
        }

        // Display 
        double norm = _sample_size * (_sample_size - 1); 
        double A_norm = static_cast<double>(A); 
        A_norm = A_norm / _sample_num / norm;
        double B_norm = static_cast<double>(B); 
        B_norm = B_norm / _sample_num / norm; 
        if (this->_output_level > 0) 
        {
            for (unsigned i = 0; i < _sample_num; ++i) 
            {
                std::cout << "[INFO] " << "A: " << ABs[2 * i] / norm 
                    << ", B: " << ABs[2 * i] / norm << std::endl; 
            }
        }
        std::cout << "A (norm): " << A_norm << ", B (norm): " 
            << B_norm << std::endl; 

        double entropy = ComputeSampen(A, B, n - K, K, this->_output_level); 
        return entropy; 
    }
private: 
    unsigned _sample_num; 
    unsigned _sample_size; 
    uniform_int_generator::random_type _rtype; 
    bool _random; 
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
    const vector<KDPoint<T, K> > points = GetKDPoints<T, K>(data_, 1);

    vector<KDPoint<T, K> > sorted_points(points);
    // The mapping p, from rank to original index 
    vector<unsigned> rank2index(n);
    for (size_t i = 0; i < n; i++) rank2index.at(i) = i;
    std::sort(rank2index.begin(), rank2index.end(),
              [&points] (unsigned i1, unsigned i2) 
              { return (points[i1] < points[i2]); });
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
        std::cout << "[INFO] The number of nodes (K = " << K << "): ";
        std::cout << tree.num_nodes() << std::endl;;
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

vector<unsigned> MergeRepeatedIndices(
    typename vector<unsigned>::const_iterator first, 
    typename vector<unsigned>::const_iterator last)
{
    size_t length = last - first; 
    vector<unsigned> counts(length, 0); 
    unsigned i = 0, k = 0; 
    while (i < length) 
    {
        unsigned count = 0; 
        while (k < length && *(first + k) == *(first + i)) 
        {
            ++k; 
            ++count; 
        }
        counts[i] = count;
        i = k; 
    }
    return counts; 
}

template<typename T, unsigned K>
long long MatchedPairsCalculatorSample<T, K>::ComputeA(
    typename vector<T>::const_iterator first,
    typename vector<T>::const_iterator last, 
    const vector<unsigned> &indices, T r)
{
    const size_t n = last - first; 
    const size_t sample_size = indices.size(); 
    vector<T> data_(first, last); 
    // Add K - 1 auxiliary points. 
    T minimum = *std::min_element(first, last); 
    // TODO: !! For debug here. 
    for (size_t i = 0; i < K - 1; i++) data_.push_back(minimum); 
    // Construct Points and merge repeated points. 
    vector<KDPoint<T, K> > points = GetKDPoints<T, K>(data_, 0); 
    for (unsigned i = 0; i < sample_size; ++i) 
    {
        unsigned index = indices[i]; 
        points[index].set_count(points[index].count() + 1); 
    }
    vector<KDPoint<T, K> > sorted_points(points.cbegin(), points.cend());
    // The mapping p, from rank to original index 
    vector<unsigned> rank2index(n); 
    for (size_t i = 0; i < n; i++) rank2index.at(i) = i; 
    std::sort(rank2index.begin(), rank2index.end(),
              [&points] (unsigned i1, unsigned i2) 
              { return (points[i1] < points[i2]); }); 
    for (size_t i = 0; i < n; i++) sorted_points[i] = points[rank2index[i]];

    MergeRepeatedPoints(sorted_points, rank2index);

    const Bounds bounds = GetRankBounds(sorted_points, r); 
    const vector<KDPoint<unsigned, K - 1> > points_grid = Map2Grid(
        sorted_points, rank2index); 

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
        long long count_repeated =\
            static_cast<long long>(points_count[i].count());
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
        std::cout << "[INFO] The number of nodes (K = " << K << "): ";
        std::cout << tree.num_nodes() << std::endl;;
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
              [&points](unsigned i1, unsigned i2) 
              { return (points[i1] < points[i2]); });
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
