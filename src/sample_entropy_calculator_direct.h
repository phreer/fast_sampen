#ifndef __SAMPLE_ENTROPY_CALCULATOR_DIRECT__
#define __SAMPLE_ENTROPY_CALCULATOR_DIRECT__

#include <string.h>
#include <vector>

#include "sample_entropy_calculator.h"
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
inline double SampleEntropyCalculatorDirect<T, K>::ComputeSampleEntropy(typename vector<T>::const_iterator first, typename vector<T>::const_iterator last, T r)
{
    long long a = 0;
    long long b = 0;

    vector<KDPoint<T, K + 1> > points = GetKDPoints<T, K + 1>(vector<T>(first, last));
    const unsigned n = points.size();
    for (unsigned i = 0; i < n; ++i)
    {
        for (unsigned j = i + 1; j < n; ++j)
        {
            if (points[i].Within(points[j], r, K))
            {
                ++b;
                T diff = points[j][K] - points[i][K];
                if (-r <= diff && diff <= r) a++;
            }
        }
    }
    return ComputeSampen(
        static_cast<double>(a),
        static_cast<double>(b),
        n - K, K, this->_output_level);
}

template<typename T, unsigned K>
inline double SampleEntropyCalculatorFastDirect<T, K>::ComputeSampleEntropy(typename vector<T>::const_iterator first, typename vector<T>::const_iterator last, T r)
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
    double result = ComputeSampen(
        static_cast<double>(ab[0]),
        static_cast<double>(ab[1]),
        n - K, r, this->_output_level);
    return result;
}

} // namespace kdtree_mddc


#endif // !__SAMPLE_ENTROPY_CALCULATOR_DIRECT__