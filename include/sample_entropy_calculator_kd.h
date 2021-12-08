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
class MatchedPairsCalculatorSampling
{
public: 
    MatchedPairsCalculatorSampling(OutputLevel output_level) 
        : _output_level(output_level) {} 
    vector<long long> ComputeA(typename vector<T>::const_iterator first, 
                               typename vector<T>::const_iterator last, 
                               unsigned sample_num, 
                               const vector<unsigned> &indices, 
                               T r); 
private: 
    OutputLevel _output_level; 
}; 


template<typename T, unsigned K> 
class MatchedPairsCalculatorSampling2
{
public: 
    MatchedPairsCalculatorSampling2(OutputLevel output_level) 
        : _output_level(output_level) {} 
    long long ComputeA(typename vector<T>::const_iterator first, 
                       typename vector<T>::const_iterator last, 
                       T r,
                       vector<unsigned>& indices); 
private: 
    OutputLevel _output_level; 
}; 


template<typename T, unsigned K>
class SampleEntropyCalculatorMao : public SampleEntropyCalculator<T, K>
{
public:
    using SampleEntropyCalculator<T, K>::SampleEntropyCalculator;
    std::string get_result_str() override
    {
        std::stringstream ss;
        ss << this->SampleEntropyCalculator<T, K>::get_result_str();
        ss << "----------------------------------------"
            << "----------------------------------------\n";
        return ss.str();
    }
protected: 
    void _ComputeSampleEntropy() override 
    {
        if (_n <= K)
        {
            std::cerr << "Data length is too short (n = " << _n;
            std::cerr << ", K = " << K << ")" << std::endl;
            exit(-1); 
        }

        MatchedPairsCalculatorMao<T, K> b_cal(this->_output_level);
        MatchedPairsCalculatorMao<T, K + 1> a_cal(this->_output_level);
        _b = b_cal.ComputeA(_data.cbegin(), _data.cend() - 1, _r);
        _a = a_cal.ComputeA(_data.cbegin(), _data.cend(), _r);
    }
    std::string _Method() const override 
    { return std::string("kd tree (Mao)"); }
    using SampleEntropyCalculator<T, K>::_data; 
    using SampleEntropyCalculator<T, K>::_r; 
    using SampleEntropyCalculator<T, K>::_n; 
    using SampleEntropyCalculator<T, K>::_computed; 
    using SampleEntropyCalculator<T, K>::_a; 
    using SampleEntropyCalculator<T, K>::_b; 
    using SampleEntropyCalculator<T, K>::_output_level; 
    using SampleEntropyCalculator<T, K>::_elapsed_seconds; 
};


template<typename T, unsigned K>
class SampleEntropyCalculatorSamplingMao : public SampleEntropyCalculatorSampling<T, K>
{
public:
    using SampleEntropyCalculatorSampling<T, K>::SampleEntropyCalculatorSampling;
    SampleEntropyCalculatorSamplingMao(
        typename vector<T>::const_iterator first, 
        typename vector<T>::const_iterator last, 
        T r, 
        unsigned sample_size, 
        unsigned sample_num, 
        double real_entropy, 
        double real_a_norm, 
        double real_b_norm, 
        RandomType rtype, 
        bool random_, 
        OutputLevel output_level) : 
        SampleEntropyCalculatorSampling<T, K>(
            first, last, r, sample_size, sample_num, real_entropy, 
            real_a_norm, real_b_norm, output_level), 
        _rtype(rtype), _random(random_) 
    {
        if (sample_num != 1) {
            std::cerr << "Only support the parameter sample_num == 1.\n" << std::endl;
            exit(-1);
        }
    }
    std::string get_result_str() override
    {
        std::stringstream ss;
        ss << this->SampleEntropyCalculatorSampling<T, K>::get_result_str();
        ss << "----------------------------------------"
            << "----------------------------------------\n";
        return ss.str();
    }
protected: 
    void _ComputeSampleEntropy() override 
    {
        if (_n <= K)
        {
            std::cerr << "Data length is too short (n = " << _n;
            std::cerr << ", K = " << K << ")" << std::endl;
            exit(-1); 
        }
        RandomIndicesSamplerWR sampler(_n - K, _sample_size, 1, _rtype, _random);
        std::vector<unsigned> sample_indices = sampler.GetSampleArrays()[0];
        MatchedPairsCalculatorSampling2<T, K> b_cal(this->_output_level);
        MatchedPairsCalculatorSampling2<T, K + 1> a_cal(this->_output_level);
        _b = b_cal.ComputeA(_data.cbegin(), _data.cend() - 1, _r, sample_indices);
        _a = a_cal.ComputeA(_data.cbegin(), _data.cend(), _r, sample_indices);
        _a_vec = std::vector<long long>(1, _a);
        _b_vec = std::vector<long long>(1, _b);
    }
    std::string _Method() const override 
    { return std::string("kd tree (Mao) sampling"); }
    using SampleEntropyCalculatorSampling<T, K>::_data; 
    using SampleEntropyCalculatorSampling<T, K>::_r; 
    using SampleEntropyCalculatorSampling<T, K>::_n; 
    using SampleEntropyCalculatorSampling<T, K>::_computed; 
    using SampleEntropyCalculatorSampling<T, K>::_a; 
    using SampleEntropyCalculatorSampling<T, K>::_b; 
    using SampleEntropyCalculatorSampling<T, K>::_output_level; 
    using SampleEntropyCalculatorSampling<T, K>::_elapsed_seconds; 
    using SampleEntropyCalculatorSampling<T, K>::_sample_size;
    using SampleEntropyCalculatorSampling<T, K>::_sample_num;
    using SampleEntropyCalculatorSampling<T, K>::_a_vec;
    using SampleEntropyCalculatorSampling<T, K>::_b_vec;
    RandomType _rtype;
    bool _random;
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


template<typename T, unsigned K> 
class ABCalculatorSamplingLiu 
{
public: 
    ABCalculatorSamplingLiu(OutputLevel output_level):
        _output_level(output_level) {} 
    vector<long long> ComputeAB(typename vector<T>::const_iterator first, 
                                typename vector<T>::const_iterator last,
                                unsigned sample_num, 
                                const vector<unsigned> &indices, 
                                T r);
private:
    OutputLevel _output_level;
}; 

/**
 * @brief This class calculates matched pairs like the class
 * MatchedPairCalculatorMao, except that only one kd tree is used.
 */
template<typename T, unsigned K>
class SampleEntropyCalculatorLiu : public SampleEntropyCalculator<T, K>
{
public:
    using SampleEntropyCalculator<T, K>::SampleEntropyCalculator;
    std::string get_result_str() override
    {
        std::stringstream ss;
        ss << this->SampleEntropyCalculator<T, K>::get_result_str();
        ss << "----------------------------------------"
            << "----------------------------------------\n";
        return ss.str();
    }
protected: 
    void _ComputeSampleEntropy() override 
    {
        if (_n <= K)
        {
            std::cerr << "Data length is too short (n = " << _n;
            std::cerr << ", K = " << K << ")" << std::endl;
            exit(-1); 
        }
        ABCalculatorLiu<T, K> abc(this->_output_level);
        vector<long long> result = abc.ComputeAB(
            _data.cbegin(), _data.cend(), _r);
        _a = result[0]; 
        _b = result[1]; 
    }
    std::string _Method() const override 
    { return std::string("kd tree (Liu)"); }

    using SampleEntropyCalculator<T, K>::_data; 
    using SampleEntropyCalculator<T, K>::_r; 
    using SampleEntropyCalculator<T, K>::_n; 
    using SampleEntropyCalculator<T, K>::_computed; 
    using SampleEntropyCalculator<T, K>::_a; 
    using SampleEntropyCalculator<T, K>::_b; 
    using SampleEntropyCalculator<T, K>::_output_level; 
    using SampleEntropyCalculator<T, K>::_elapsed_seconds; 
};


template<typename T, unsigned K>
class SampleEntropyCalculatorSamplingKDTree : 
    public SampleEntropyCalculatorSampling<T, K> 
{
public:
    SampleEntropyCalculatorSamplingKDTree(
        typename vector<T>::const_iterator first, 
        typename vector<T>::const_iterator last, 
        T r, 
        unsigned sample_size, 
        unsigned sample_num, 
        double real_entropy, 
        double real_a_norm, 
        double real_b_norm, 
        RandomType rtype, 
        bool random_, 
        OutputLevel output_level) : 
        SampleEntropyCalculatorSampling<T, K>(
            first, last, r, sample_size, sample_num, real_entropy, 
            real_a_norm, real_b_norm, output_level), 
        _rtype(rtype), _random(random_) 
    {}
    std::string get_result_str() override 
    {
        std::stringstream ss; 
        ss << this->SampleEntropyCalculatorSampling<T, K>::get_result_str(); 
        if (_output_level) 
        {
            ss.precision(4); 
            ss << std::scientific; 
            ss << "[INFO] " << "random_type: " 
                << random_type_names[_rtype] << "\n" 
                << "[INFO] " << "random: " << _random << "\n"; 
        }
        ss << "----------------------------------------"
            << "----------------------------------------\n";
        return ss.str(); 
    }

protected:
    void _ComputeSampleEntropy() override 
    {
        if (_n <= K) 
        {
            std::cerr << "Data length is too short (n = " << _n; 
            std::cerr << ", K = " << K << ")" << std::endl; 
            exit(-1); 
        }
        RandomIndicesSampler uig(0, _n - 2, _rtype, _random); 
        vector<unsigned> indices(_sample_size * _sample_num); 
        for (unsigned i = 0; i < _sample_num; ++i) 
        {
            for (unsigned j = 0; j < _sample_size; ++j) 
            {
                indices[i * _sample_size + j] = uig.get(); 
            }
        }
        
        MatchedPairsCalculatorSampling<T, K> b_cal(this->_output_level);
        MatchedPairsCalculatorSampling<T, K + 1> a_cal(this->_output_level);
        _a_vec = a_cal.ComputeA(
            _data.cbegin(), _data.cend(), _sample_num, indices, _r);
        _b_vec = b_cal.ComputeA(
            _data.cbegin(), _data.cend(), _sample_num, indices, _r);
    }

    std::string _Method() const override 
    { 
        return std::string("sampling kd tree (Mao, ") 
            + random_type_names[_rtype] + std::string(")"); 
    }

    using SampleEntropyCalculatorSampling<T, K>::_data; 
    using SampleEntropyCalculatorSampling<T, K>::_r; 
    using SampleEntropyCalculatorSampling<T, K>::_n; 
    using SampleEntropyCalculatorSampling<T, K>::_computed; 
    using SampleEntropyCalculatorSampling<T, K>::_a; 
    using SampleEntropyCalculatorSampling<T, K>::_b; 
    using SampleEntropyCalculatorSampling<T, K>::_output_level; 
    using SampleEntropyCalculatorSampling<T, K>::_elapsed_seconds; 
    using SampleEntropyCalculatorSampling<T, K>::_sample_size; 
    using SampleEntropyCalculatorSampling<T, K>::_sample_num; 
    using SampleEntropyCalculatorSampling<T, K>::_a_vec; 
    using SampleEntropyCalculatorSampling<T, K>::_b_vec; 
    RandomType _rtype;
    bool _random; 
};


template<typename T, unsigned K>
class SampleEntropyCalculatorSamplingLiu :
    public SampleEntropyCalculatorSampling<T, K>
{
public:
    SampleEntropyCalculatorSamplingLiu(
        typename vector<T>::const_iterator first,
        typename vector<T>::const_iterator last,
        T r,
        unsigned sample_size,
        unsigned sample_num,
        double real_entropy,
        double real_a_norm,
        double real_b_norm,
        RandomType rtype,
        bool random_,
        OutputLevel output_level) :
        SampleEntropyCalculatorSampling<T, K>(
            first, last, r, sample_size, sample_num, real_entropy,
            real_a_norm, real_b_norm, output_level),
        _rtype(rtype), _random(random_)
    {}
    std::string get_result_str() override
    {
        std::stringstream ss;
        ss << this->SampleEntropyCalculatorSampling<T, K>::get_result_str();
        if (_output_level)
        {
            ss.precision(4);
            ss << std::scientific;
            ss << "[INFO] " << "random_type: "
                << random_type_names[_rtype] << "\n"
                << "[INFO] " << "random: " << _random << "\n";
        }
        ss << "----------------------------------------"
            << "----------------------------------------\n";
        return ss.str();
    }

protected:
    void _ComputeSampleEntropy() override
    {
        if (_n <= K)
        {
            std::cerr << "Data length is too short (n = " << _n;
            std::cerr << ", K = " << K << ")" << std::endl;
            exit(-1);
        }
        RandomIndicesSampler uig(0, _n - 2, _rtype, _random);
        vector<unsigned> indices(_sample_size * _sample_num);
        for (unsigned i = 0; i < _sample_num; ++i)
        {
            for (unsigned j = 0; j < _sample_size; ++j)
            {
                indices[i * _sample_size + j] = uig.get();
            }
        }
        
        ABCalculatorSamplingLiu<T, K> ab_cal(_output_level);
        vector<long long> results = ab_cal.ComputeAB(
            _data.cbegin(), _data.cend(), _sample_num, indices, _r);
        
        _a_vec = vector<long long>(_sample_num);
        _b_vec = vector<long long>(_sample_num);
        for (unsigned i = 0; i < _sample_num; ++i)
        {
            _a_vec[i] = results[2 * i];
            _b_vec[i] = results[2 * i + 1];
        }
    }
    std::string _Method() const override
    {
        return std::string("sampling kd tree (liu, ")
            + random_type_names[_rtype] + std::string(")");
    }

    using SampleEntropyCalculatorSampling<T, K>::_data;
    using SampleEntropyCalculatorSampling<T, K>::_r;
    using SampleEntropyCalculatorSampling<T, K>::_n;
    using SampleEntropyCalculatorSampling<T, K>::_computed;
    using SampleEntropyCalculatorSampling<T, K>::_a;
    using SampleEntropyCalculatorSampling<T, K>::_b;
    using SampleEntropyCalculatorSampling<T, K>::_output_level;
    using SampleEntropyCalculatorSampling<T, K>::_elapsed_seconds;
    using SampleEntropyCalculatorSampling<T, K>::_sample_size;
    using SampleEntropyCalculatorSampling<T, K>::_sample_num;
    using SampleEntropyCalculatorSampling<T, K>::_a_vec;
    using SampleEntropyCalculatorSampling<T, K>::_b_vec;
    RandomType _rtype;
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
    const vector<KDPoint<T, K> > points = GetKDPoints<T, K>(
        data_.cbegin(), data_.cend(), 1);

    vector<KDPoint<T, K> > sorted_points(points);
    // The mapping p, from rank to original index 
    vector<unsigned> rank2index(n);
    for (size_t i = 0; i < n; i++) rank2index.at(i) = i;

    Timer timer; 
    std::sort(rank2index.begin(), rank2index.end(),
              [&points] (unsigned i1, unsigned i2) 
              { return (points[i1] < points[i2]); });
    timer.StopTimer(); 
    if (_output_level == Debug) 
    {
        std::cout << "[INFO] Time consumed in presorting: " 
            << timer.ElapsedSeconds() << "s\n";
    }

    for (size_t i = 0; i < n; i++) sorted_points[i] = points[rank2index[i]];

    MergeRepeatedPoints(sorted_points, rank2index);

    const Bounds bounds = GetRankBounds(sorted_points, r);
    const vector<KDPoint<unsigned, K - 1> > points_grid = 
        Map2Grid(sorted_points, rank2index);

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
    
    timer.SetStartingPointNow();
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
    timer.StopTimer(); 
    
    if (_output_level == Debug)
    {
        std::cout << "[INFO] Time consumed in range counting: "
            << timer.ElapsedSeconds() << " seconds\n";
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
long long MatchedPairsCalculatorSampling2<T, K>::ComputeA(
    typename vector<T>::const_iterator first,
    typename vector<T>::const_iterator last,
    T r, std::vector<unsigned>& sample_indices)
{
    const size_t n = last - first;
    assert(n - K + 1 >= sample_indices.size());
    vector<T> data_(first, last);
    // Add K - 1 auxiliary points. 
    T minimum = *std::min_element(first, last);
    for (size_t i = 0; i < K - 1; ++i) data_.push_back(minimum);

    // Construct points in k-dimensional space and merge repeated points. 
    vector<KDPoint<T, K> > points = GetKDPoints<T, K>(
        data_.cbegin(), data_.cend(), 1);
    for (size_t i = points.size() - K + 1; i < points.size(); ++i) {
        points[i].set_count(0);
    }

    vector<KDPoint<T, K> > sorted_points(points);
    // The mapping p, from rank to original index 
    vector<unsigned> rank2index(n);
    for (size_t i = 0; i < n; i++) rank2index.at(i) = i;
    Timer timer; 
    std::sort(rank2index.begin(), rank2index.end(),
              [&points] (unsigned i1, unsigned i2) 
              { return (points[i1] < points[i2]); });
    timer.StopTimer(); 
    if (_output_level == Debug) 
    {
        std::cout << "[INFO] Time consumed in presorting: " 
            << timer.ElapsedSeconds() << "s\n";
    }

    for (size_t i = 0; i < n; i++) sorted_points[i] = points[rank2index[i]];

    const Bounds bounds = GetRankBounds(sorted_points, r);
    // Map values of each coordinate to the rank given by sorting.
    // Since the value at first dimension equal to the index of that point
    // in the sorted array, we can reduce the dimension of the points by 1.
    const vector<KDPoint<unsigned, K - 1> > points_grid = 
        Map2Grid(sorted_points, rank2index);

    // Construct kd tree.
    vector<KDPoint<unsigned, K - 1> > points_count;
    vector<unsigned> points_count_indices;
    // Only use points with positive count to construct the kd tree.
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
    assert(n_count == n - K + 1);
    int i_sample_indices = 0;
    std::sort(sample_indices.begin(), sample_indices.end());
    timer.SetStartingPointNow();
    for (unsigned i = 0; i < n_count; i++)
    {
        if (sample_indices[i_sample_indices] != i) {
            continue;
        }         
        ++i_sample_indices;

        const unsigned rank1 = points_count_indices[i];
        unsigned upperbound = bounds.upper_bounds[rank1];

        if (upperbound < points_count_indices[i + 1]) continue;
        // Close nodes whose value of the first dimension are outside bounds.
        if (i_sample_indices > 0) {
            for (unsigned k = sample_indices[i_sample_indices - 1] + 1; k <= i; ++k) {
                tree.Close(k);
            }
        }

        // Update tree. 
        if (upperbound_prev < rank1) upperbound_prev = rank1;
        unsigned j = i + 1;
        while (j < n_count && points_count_indices[j] <= upperbound_prev) ++j;
        while (j < n_count && points_count_indices[j] <= upperbound) {
            tree.UpdateCount(j, points_count[j].count());
            ++num_opened;
            ++j;
        }

        const Range<unsigned, K - 1> range = GetHyperCube(points_count[i], bounds);
        result += tree.CountRange(range, num_nodes);
        ++num_countrange_called;
        upperbound_prev = upperbound;
    }
    timer.StopTimer(); 

    std::cout << "i_sample_indices: " << i_sample_indices << std::endl;
    assert(i_sample_indices == sample_indices.size());

    if (_output_level == Debug)
    {
        std::cout << "[INFO] Time consumed in range counting: "
            << timer.ElapsedSeconds() << " seconds\n";
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
vector<long long> MatchedPairsCalculatorSampling<T, K>::ComputeA(
    typename vector<T>::const_iterator first,
    typename vector<T>::const_iterator last, 
    unsigned sample_num, 
    const vector<unsigned> &indices, T r)
{
    const size_t n = last - first; 
    const size_t sample_size = indices.size() / sample_num; 
    vector<T> data_(first, last); 
    // Add K - 1 auxiliary points. 
    T minimum = *std::min_element(first, last); 
    // TODO: !! For debug here. 
    for (size_t i = 0; i < K - 1; i++) data_.push_back(minimum); 
    // Construct Points and merge repeated points. 
    vector<KDPoint<T, K> > points = GetKDPoints<T, K>(data_.cbegin(),
                                                      data_.cend(), 0); 
    // The mapping p, from rank to original index 
    vector<unsigned> rank2index(n); 
    for (size_t i = 0; i < n; i++) rank2index.at(i) = i; 

    Timer timer; 
    std::sort(rank2index.begin(), rank2index.end(),
              [&points] (unsigned i1, unsigned i2) 
              { return (points[i1] < points[i2]); }); 
    timer.StopTimer(); 
    if (_output_level == Debug) 
    {
        std::cout << "[INFO] Time consumed in presorting: " 
            << timer.ElapsedSeconds() << "s\n";
    }

    vector<KDPoint<T, K> > sorted_points(points.cbegin(), points.cend());
    for (size_t i = 0; i < n; i++) 
        sorted_points[i] = points[rank2index[i]];

    // MergeRepeatedPoints(sorted_points, rank2index);

    const Bounds bounds = GetRankBounds(sorted_points, r); 
    vector<KDPoint<unsigned, K - 1> > points_grid = Map2Grid(
        sorted_points, rank2index, false);

    vector<vector<unsigned> > sample_groups(n); 
    for (unsigned i = 0; i < sample_num * sample_size; ++i) 
    {
        unsigned index = indices[i];
        points_grid[index].increase_count(1);
        sample_groups[index].push_back(i / sample_size); 
    }

    // Construct kd tree.
    vector<KDPoint<unsigned, K - 1> > points_count;
    vector<unsigned> points_count_indices;
    vector<vector<unsigned> > indices_count(sample_num); 
    for (unsigned i = 0; i < n; i++)
    {
        if (points_grid[i].count())
        {
            for (unsigned j = 0; j < sample_groups[i].size(); ++j) 
            {
                indices_count[sample_groups[i][j]].
                    push_back(points_count.size());
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
    for (unsigned i = 0; i < n_count; ++i) 
    {
        points_count[i].set_count(0); 
    }
    vector<long long> results(sample_num); 
    for (unsigned k = 0; k < sample_num; ++k) 
    {
        for (unsigned i = 0; i < indices_count[k].size(); ++i) 
        {
            points_count[indices_count[k][i]].increase_count(1); 
        }
        long long result = 0;
        for (unsigned i = 0; i < n_count - 1; ++i)
        {
            if (points_count[i].count() == 0) 
                continue; 

            const unsigned rank1 = points_count_indices[i];
            unsigned upperbound = bounds.upper_bounds[rank1];
            long long count_repeated =
                static_cast<long long>(points_count[i].count());
            result += (count_repeated - 1) * count_repeated / 2;

            if (upperbound < points_count_indices[i + 1])
                continue;
            // Close current node.
            tree.Close(i);

            // Update tree.
            if (upperbound_prev < rank1)
                upperbound_prev = rank1;
            unsigned j = i + 1;
            while (j < n_count && points_count_indices[j] <= upperbound_prev)
                ++j;
            while (j < n_count && points_count_indices[j] <= upperbound)
            {
                if (points_count[j].count()) 
                {
                    tree.UpdateCount(j, points_count[j].count());
                    ++num_opened;
                }
                ++j;
            }

            const Range<unsigned, K - 1> range = GetHyperCube(points_count[i], bounds);
            result += tree.CountRange(range, num_nodes) * count_repeated;
            ++num_countrange_called;
            upperbound_prev = upperbound;
        }
        tree.Close(n_count - 1); 
        results[k] = result; 

        // Close points for current sample group. 
        for (unsigned i = 0; i < indices_count[k].size(); ++i) 
        {
            points_count[indices_count[k][i]].set_count(0); 
        }
    }
    timer.StopTimer();
    
    if (_output_level == Debug)
    {
        std::cout << "[INFO] Time consumed in range counting: "
            << timer.ElapsedSeconds() << " seconds\n";
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
    if (_output_level >= 2) 
    {
        std::cout << "[DEBUG] The numbers of the matched pairs: \n";
        for (unsigned i = 0; i < sample_num; ++i) 
        {
            if (i) std::cout << ", " << results[i]; 
            else std::cout << results[i]; 
        }
        std::cout << std::endl; 
        for (unsigned i = 0; i < sample_num; ++i) 
        {
            std::cout << "[DEBUG] Sample indices (count): \n"; 
            for (unsigned j = 0; j < sample_size; ++j) 
            {
                std::cout << indices_count[i][j] << ", "; 
            }
            std::cout << std::endl; 
        }
    }
    return results;
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
    const vector<KDPoint<T, K + 1> > points = GetKDPoints<T, K + 1>(
        data_.cbegin(), data_.cend());

    vector<KDPoint<T, K + 1> > sorted_points(points);
    // The mapping p, from rank to original index 
    vector<unsigned> rank2index(n);
    for (size_t i = 0; i < n; i++) rank2index.at(i) = i;

    Timer timer; 
    std::sort(rank2index.begin(), rank2index.end(),
              [&points](unsigned i1, unsigned i2) 
              { return (points[i1] < points[i2]); });
    for (size_t i = 0; i < n; i++) sorted_points[i] = points[rank2index[i]];
    timer.StopTimer(); 
    if (_output_level == Debug)
    {
        std::cout << "[INFO] Time consumed in presorting: " 
            << timer.ElapsedSeconds() << " seconds\n";
    }

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
    
    timer.SetStartingPointNow();
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
    timer.StopTimer();

    if (_output_level == Debug)
    {
        std::cout << "[INFO] Time consumed in range counting: "
            << timer.ElapsedSeconds() << " seconds\n";
        std::cout << "[INFO] The number of nodes (K = " << K << "): ";
        std::cout << tree.num_nodes() << std::endl;;
        std::cout << "[INFO] The number of leaf nodes (K = " << K << "): ";
        std::cout << n_count << std::endl;
        std::cout << "[INFO] The number of calls for CountRange(): ";
        std::cout << num_countrange_called << std::endl;
        std::cout << "[INFO] The number times to open node: ";
        std::cout << num_opened << std::endl;
        std::cout << "[INFO] The number of nodes visited (K = " << K << "): ";
        std::cout << num_nodes << std::endl;
    }

    return result;
}


template<typename T, unsigned K>
vector<long long> 
ABCalculatorSamplingLiu<T, K>::ComputeAB(
    typename vector<T>::const_iterator first, 
    typename vector<T>::const_iterator last, 
    unsigned sample_num, 
    const vector<unsigned> &indices, 
    T r)
{
    const unsigned n = last - first;
    vector<T> data_(first, last);
    const unsigned sample_size = indices.size() / sample_num;
    // Add K auxiliary points. 
    T minimum = *std::min_element(first, last);
    for (size_t i = 0; i < K; i++) data_.push_back(minimum);
    // Construct Points and merge repeated points. 
    const vector<KDPoint<T, K + 1> > points = GetKDPoints<T, K + 1>(
        data_.cbegin(), data_.cend(), 0);

    vector<KDPoint<T, K + 1> > sorted_points(points);
    // The mapping p, from rank to original index 
    vector<unsigned> rank2index(n);
    for (size_t i = 0; i < n; i++) rank2index.at(i) = i;

    Timer timer; 
    std::sort(rank2index.begin(), rank2index.end(),
              [&points](unsigned i1, unsigned i2) 
              { return (points[i1] < points[i2]); });
    for (size_t i = 0; i < n; i++) sorted_points[i] = points[rank2index[i]];
    timer.StopTimer(); 
    if (_output_level == Debug)
    {
        std::cout << "[INFO] Time consumed in presorting: " 
            << timer.ElapsedSeconds() << "s\n";
    }

    const Bounds bounds = GetRankBounds(sorted_points, r);
    vector<KDPoint<unsigned, K> > points_grid = Map2Grid(
        sorted_points, rank2index, false);

    vector<vector<unsigned> > sample_groups(n); 
    for (unsigned i = 0; i < sample_num * sample_size; ++i)
    {
        unsigned index = indices[i];
        points_grid[index].increase_count(1);
        sample_groups[index].push_back(i / sample_size); 
    }

    // Construct kd tree.
    vector<KDPoint<unsigned, K> > points_count;
    vector<unsigned> points_count_indices;
    vector<vector<unsigned> > indices_count(sample_num); 
    for (unsigned i = 0; i < n; i++)
    {
        if (points_grid[i].count())
        {
            for (unsigned j = 0; j < sample_groups[i].size(); ++j) 
            {
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
    for (unsigned i = 0; i < n_count; ++i) 
    {
        points_count[i].set_count(0); 
    }
    for (unsigned k = 0; k < sample_num; ++k) 
    {
        for (unsigned i = 0; i < indices_count[k].size(); ++i) 
        {
            points_count[indices_count[k][i]].increase_count(1); 
        }
        long long result_a = 0; 
        long long result_b = 0; 
        for (unsigned i = 0; i < n_count - 1; ++i)
        {
            if (points_count[i].count() == 0) continue; 

            const unsigned rank1 = points_count_indices[i];
            unsigned upperbound = bounds.upper_bounds[rank1];
            long long count_repeated = 
                static_cast<long long>(points_count[i].count());
            result_a += (count_repeated - 1) * count_repeated / 2; 
            result_b += (count_repeated - 1) * count_repeated / 2; 

            if (upperbound < points_count_indices[i + 1]) continue;
            // Close current node. 
            tree.Close(i);

            // Update tree. 
            if (upperbound_prev < rank1) upperbound_prev = rank1;
            unsigned j = i + 1;
            while (j < n_count && points_count_indices[j] <= upperbound_prev) 
                ++j;
            while (j < n_count && points_count_indices[j] <= upperbound)
            {
                if (points_count[j].count())
                {
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
        results[2 * k + 1]= result_b;
        tree.Close(n_count - 1); 

        for (unsigned i = 0; i < indices_count[k].size(); ++i) 
        {
            points_count[indices_count[k][i]].set_count(0); 
        }
    }
    timer.StopTimer();
    
    if (_output_level == Debug)
    {
        std::cout << "[INFO] Time consumed in range counting: "
            << timer.ElapsedSeconds() << " seconds\n";
        std::cout << "[INFO] The number of nodes (K = " << K << "): ";
        std::cout << tree.num_nodes() << std::endl;;
        std::cout << "[INFO] The number of leaf nodes (K = " << K << "): ";
        std::cout << n_count << std::endl;
        std::cout << "[INFO] The number of calls for CountRange(): ";
        std::cout << num_countrange_called << std::endl;
        std::cout << "[INFO] The number times to open node: ";
        std::cout << num_opened << std::endl;
        std::cout << "[INFO] The number of nodes visited (K = " << K << "): ";
        std::cout << num_nodes << std::endl;
    }

    return results; 
}

} // namespace kdtree_mddc

#endif // !__SAMPLE_ENTROPY_CALCULATOR_KD__
