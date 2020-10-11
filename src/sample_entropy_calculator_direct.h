#ifndef __SAMPLE_ENTROPY_CALCULATOR_DIRECT__
#define __SAMPLE_ENTROPY_CALCULATOR_DIRECT__

#include <string.h>
#include <vector>
#include <numeric> 

#include "sample_entropy_calculator.h"
#include "random_sampler.h"

namespace kdtree_mddc
{

template<typename T, unsigned K>
vector<long long> _ComputeABFastDirect(const T *y, unsigned n, T r);

template<typename T, unsigned K>
class SampleEntropyCalculatorDirect : public SampleEntropyCalculator<T, K>
{
public:
    SampleEntropyCalculatorDirect(typename vector<T>::const_iterator first,
                                  typename vector<T>::const_iterator last, 
                                  T r, OutputLevel output_level) : 
        SampleEntropyCalculator<T, K>(first, last, r, output_level) 
    {}
    std::string get_result_str() override
    {
        std::stringstream ss;
        ss << this->SampleEntropyCalculator<T, K>::get_result_str();
        ss << "----------------------------------------"
            << "----------------------------------------\n";
        return ss.str();
    }
protected: 
    void _ComputeSampleEntropy() override; 
    std::string _Method() const override 
    { return std::string("plain direct"); }
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
class SampleEntropyCalculatorFastDirect : public SampleEntropyCalculator<T, K>
{
public:
    SampleEntropyCalculatorFastDirect(typename vector<T>::const_iterator first,
                                      typename vector<T>::const_iterator last, 
                                      T r, OutputLevel output_level) : 
        SampleEntropyCalculator<T, K>(first, last, r, output_level) 
    {}
    std::string get_result_str() override
    {
        std::stringstream ss;
        ss << this->SampleEntropyCalculator<T, K>::get_result_str();
        ss << "----------------------------------------"
            << "----------------------------------------\n";
        return ss.str();
    }
protected: 
    void _ComputeSampleEntropy() override;
    std::string _Method() const override 
    { return std::string("fast direct"); }
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
class SampleEntropyCalculatorSamplingDirect : 
    public SampleEntropyCalculatorSampling<T, K>
{
public:
    SampleEntropyCalculatorSamplingDirect(
        typename vector<T>::const_iterator first,
        typename vector<T>::const_iterator last, 
        T r, unsigned sample_size, unsigned sample_num,
        double real_entropy, double real_a_norm, double real_b_norm, 
        RandomType rtype, bool random_, OutputLevel output_level) : 
        SampleEntropyCalculatorSampling<T, K>(
            first, last, r, sample_size, sample_num, real_entropy, 
            real_a_norm, real_b_norm, output_level), 
        _rtype(rtype), _random(random_) 
    {}
    std::string get_result_str() override
    {
        std::stringstream ss;
        ss << this->SampleEntropyCalculatorSampling<T, K>::get_result_str();
        ss << "----------------------------------------"
            << "----------------------------------------\n";
        return ss.str();
    }
protected:
    void _ComputeSampleEntropy() override;
    std::string _Method() const override 
    { 
        return std::string("sampling direct (") 
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


/*
* Calculates an estimate of sample entropy but does NOT calculate
* the variance of the estimate
*/
template<typename T, unsigned K>
vector<long long> _ComputeABFastDirect(const T *y, unsigned n, T r)
{
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
    for (unsigned i = 0; i < n - 1; i++)
    {
        nj = n - i - 1;
        y1 = y[i];
        for (jj = 0; jj < nj; jj++)
        {
            j = jj + i + 1;
            if (((y[j] - y1) <= r) && ((y1 - y[j]) <= r))
            {
                run[jj] = lastrun[jj] + 1;
                M1 = M < run[jj] ? M : run[jj];
                for (m = 0; m < M1; m++)
                {
                    A[m]++;
                    if (j < n - 1)
                        B[m]++;
                }
            } else
                run[jj] = 0;
        }			/* for jj */
        for (j = 0; j < nj; j++)
            lastrun[j] = run[j];
    }				/* for i */

    delete[] run;
    delete[] lastrun;

    vector<long long> result(2, 0);
    result[0] = A[M - 1];
    result[1] = B[M - 2];
    return result;
}

template<typename T, unsigned K>
vector<long long> ComputeABDirect(const vector<KDPoint<T, K + 1> >&points, T r)
{
    const unsigned n = points.size();
    vector<long long> results(2); 
    long long a = 0LL, b = 0LL; 
    for (unsigned i = 0; i < n; ++i)
    {
        for (unsigned j = i + 1; j < n; ++j)
        {
            if (points[i].Within(points[j], r, K))
            {
                ++b;
                T diff = points[j][K] - points[i][K];
                if (-r <= diff && diff <= r) ++a;
            }
        }
    }
    results[0] = a; 
    results[1] = b;
    return results;
}

template<typename T, unsigned K>
void SampleEntropyCalculatorDirect<T, K>::_ComputeSampleEntropy()
{
    vector<KDPoint<T, K + 1> > points = GetKDPoints<T, K + 1>(
        SampleEntropyCalculator<T, K>::_data.cbegin(), 
        SampleEntropyCalculator<T, K>::_data.cend());
    vector<long long> ab = ComputeABDirect<T, K>(
        points, SampleEntropyCalculator<T, K>::_r);
    SampleEntropyCalculator<T, K>::_a = ab[0]; 
    SampleEntropyCalculator<T, K>::_b = ab[1]; 
}

template<typename T, unsigned K>
void SampleEntropyCalculatorFastDirect<T, K>::_ComputeSampleEntropy()
{
    if (_n <= K)
    {
        std::cerr << "Data length is too short (n = " << _n;
        std::cerr << ", K = " << K << ")" << std::endl;
        exit(-1); 
    }

    vector<long long> ab = _ComputeABFastDirect<T, K>(
        _data.data(), _data.size(), _r);
    _a = ab[0], _b = ab[1]; 
}


template<typename T, unsigned K>
void SampleEntropyCalculatorSamplingDirect<T, K>::_ComputeSampleEntropy()
{
    const vector<vector<unsigned> > indices = GetSampleIndices(
        _rtype, _n - K, _sample_size, _sample_num, _random);
    _a_vec = vector<long long>(_sample_num); 
    _b_vec = vector<long long>(_sample_num); 
    vector<vector<KDPoint<T, K + 1> > > points = GetKDPointsSample<T, K + 1>(
        _data.cbegin(), _data.cend(), indices, 1); 
    for (unsigned i = 0; i < _sample_num; ++i) 
    {
        vector<long long> ab = ComputeABDirect<T, K>(points[i], _r);
        _a_vec[i] = ab[0], _b_vec[i] = ab[1]; 
        _a += ab[0]; 
        _b += ab[1]; 
    }
}
} // namespace kdtree_mddc


#endif // !__SAMPLE_ENTROPY_CALCULATOR_DIRECT__
