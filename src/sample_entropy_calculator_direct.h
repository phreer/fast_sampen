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
class SampleEntropyCalculatorDirect : public SampleEntropyCalculator<T, K>
{
public:
    SampleEntropyCalculatorDirect(OutputLevel output_level) 
        : SampleEntropyCalculator<T, K>(output_level) {}
    double ComputeSampleEntropy(typename vector<T>::const_iterator first,
                                typename vector<T>::const_iterator last, T r);
};

template<typename T, unsigned K>
vector<long long> _ComputeABFastDirect(const T *y, unsigned n, T r);


template<typename T, unsigned K>
class SampleEntropyCalculatorFastDirect : public SampleEntropyCalculator<T, K>
{
public:
    SampleEntropyCalculatorFastDirect(OutputLevel output_level) 
        : SampleEntropyCalculator<T, K>(output_level) {}
    double ComputeSampleEntropy(typename vector<T>::const_iterator first,
                                typename vector<T>::const_iterator last, T r);
};

template<typename T, unsigned K>
class SampleEntropyCalculatorDirectSample : public SampleEntropyCalculator<T, K>
{
public:
    SampleEntropyCalculatorDirectSample(unsigned sample_size, unsigned sample_num,
                                        RandomType rtype, bool random_,
                                        OutputLevel output_level)
        : SampleEntropyCalculator<T, K>(output_level), _sample_size(sample_size),
        _sample_num(sample_num), _rtype(rtype), _random(random_) {}
    double ComputeSampleEntropy(typename vector<T>::const_iterator first,
                                typename vector<T>::const_iterator last, T r);
private:
    unsigned _sample_size;
    unsigned _sample_num;
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
inline double SampleEntropyCalculatorDirect<T, K>::ComputeSampleEntropy(
    typename vector<T>::const_iterator first, 
    typename vector<T>::const_iterator last, T r)
{
    const unsigned n = last - first;
    vector<KDPoint<T, K + 1> > points = GetKDPoints<T, K + 1>(first, last);
    vector<long long> ab = ComputeABDirect<T, K>(points, r);
    long long a = ab[0], b = ab[1]; 
    double norm = (n - K + 1) * (n - K); 
    std::cout << "A (norm): " << a / norm << ", B (norm): " << b / norm 
        << std::endl;
    return ComputeSampen(
        static_cast<double>(a),
        static_cast<double>(b),
        n - K, K, this->_output_level);
}

template<typename T, unsigned K>
inline double SampleEntropyCalculatorFastDirect<T, K>::ComputeSampleEntropy(
    typename vector<T>::const_iterator first, 
    typename vector<T>::const_iterator last, T r)
{
    const unsigned n = first - last;
    if (n <= K)
    {
        std::cerr << "Data length is too short (n = " << n;
        std::cerr << ", K = " << K << ")" << std::endl;
        return -1;
    }

    vector<T> data(first, last);
    vector<long long> ab = _ComputeABFastDirect<T, K>(data.data(), data.size(), r);
    double norm = (n - K + 1) * (n - K); 
    std::cout << "A (norm): " << ab[0] / norm << ", B (norm): " << ab[1] / norm 
        << std::endl;
    double result = ComputeSampen(
        static_cast<double>(ab[0]),
        static_cast<double>(ab[1]),
        n - K, r, this->_output_level);
    return result;
}


template<typename T, unsigned K>
double SampleEntropyCalculatorDirectSample<T, K>::ComputeSampleEntropy(
    typename vector<T>::const_iterator first,
    typename vector<T>::const_iterator last, T r)
{
    vector<long long> results(2 * _sample_num);
    const unsigned n = last - first;
    long long a = 0LL, b = 0LL;
    vector<vector<unsigned> > indices = GetSampleIndices(
        _rtype, n - K, _sample_size, _sample_num, _random);
    for (unsigned i = 0; i < _sample_num; ++i) 
    {
        vector<KDPoint<T, K + 1> > points = GetKDPointsSample<T, K + 1>(
            first, last, indices[i], 1); 
        vector<long long> ab = ComputeABDirect<T, K>(points, r);
        results[2 * i] = ab[0]; 
        results[2 * i + 1] = ab[1]; 
        a += ab[0]; 
        b += ab[1]; 
    }
    double entropy = ComputeSampen(static_cast<double>(a), 
                                   static_cast<double>(b), n - K, K, 
                                   this->_output_level); 
    return entropy; 
}
} // namespace kdtree_mddc


#endif // !__SAMPLE_ENTROPY_CALCULATOR_DIRECT__
