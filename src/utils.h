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
#include <iterator>
#include <vector>
#include <string>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <assert.h>

#include "kdpoint.h"

#ifdef _DEBUG
#define DEBUG
#endif


namespace kdtree_mddc
{
using std::vector;
using std::ifstream;
using std::cerr;
using std::string;

enum OutputLevel
{
    // Only show essential result, i.e., necessary results. 
    Silent,
    // Show helpful informations, illustrating the status of functions. 
    Info,
    // Show anything that may be useful for debugging. 
    Debug
};

struct Bounds
{
    Bounds(size_t n) :lower_bounds(n), upper_bounds(n) {}
    vector<unsigned> lower_bounds;
    vector<unsigned> upper_bounds;
};

template<typename T, unsigned K>
struct Range
{
    T lower_ranges[K];
    T upper_ranges[K];
};


vector<unsigned> GetInverseMap(const vector<unsigned> &map);

double ComputeSampen(double A, double B, unsigned N, unsigned m, OutputLevel output_level);

/**
 * @brief Read data from file. 
 * 
 * @param filename: The name of the file to read. 
 * @return Read data. 
 */
template<typename T>
vector<T> ReadData(std::string filename);

template<typename T>
double ComputeVarience(const vector<T> &data);

template<typename T>
T ComputeSum(const vector<T> &data);

/**
* @brief Merge repeated points by setting count.
*
* @param[out] points: The k-dimensional points that has been sorted.
* @note Auxiliary points will be disabled.
*/

template<typename T, unsigned K>
void MergeRepeatedPoints(vector<KDPoint<T, K> > &points,
                         const vector<unsigned> &rank2index);


template<typename T, unsigned K>
vector<KDPoint<T, K> > GetKDPoints(const vector<T> &data);

template<typename T, unsigned K>
void MergeRepeatedPoints(vector<KDPoint<T, K> > &points,
                         const vector<unsigned> &rank2index);

template<typename T, unsigned K>
void CloseAuxiliaryPoints(vector<KDPoint<T, K> > &points,
                          const vector<unsigned> &rank2index);

/*
* @brief Maps the points to grids.
*
* For each point p, p[i] is mapped to the rank of p[i] for i = 1, 2, ..., K-1
*/
template<typename T, unsigned K>
vector<KDPoint<unsigned, K - 1> > Map2Grid(
    const vector<KDPoint<T, K> > &points, const vector<unsigned> &rank2index);

/**
* @brief Get the bounds of indices such that within the bound the values
* are within the threshold r.
*
* @param data: The data to find lower bounds and upper bounds.
* @note The data must have been sorted in ascending lexicographic order.
* @param r: The threshold.
* @return The bounds of the indices.
*/
template<typename T, unsigned K>
Bounds GetRankBounds(const vector<KDPoint<T, K> > &points, T r);

/*
* @brief Given a point (in grid), get the bound.
*/
template<unsigned K>
Range<unsigned, K> GetHyperCube(const KDPoint<unsigned, K> &point, const Bounds &bounds);


class ArgumentParser
{
public:
    ArgumentParser(int argc, char *argv[]) : arg_list(argv, argv + argc) {}
    string getArg(const string &arg);
    bool isOption(const string &opt);
private:
    vector<string> arg_list;
};

//////////////////////////////////////////////////////////////////////////
// Implementation 
//////////////////////////////////////////////////////////////////////////

template<typename T>
vector<T> ReadData(std::string filename)
{
    ifstream ifs(filename);
    vector<T> result;
    if (!ifs.is_open())
    {
        cerr << "Cannot open file! (filename: " << filename << ")" << std::endl;
    } else
    {
        T x = 0;
        while (ifs >> x) result.push_back(x);
        ifs.close();
    }
    return result;
}

template<typename T, unsigned K>
vector<KDPoint<T, K> >GetKDPoints(const vector<T> &data)
{
    const size_t n = data.size();
    
    if (n <= K)
    {
        std::cerr << "GetKDPoints(): Data length is too short (n = " << n;
        std::cerr << ", K = " << K << ")" << std::endl;
        return vector<KDPoint<T, K> >();
    } else
    {
        vector<KDPoint<T, K> > points(n - K + 1);
        auto begin_i = data.cbegin();
        for (size_t i = 0; i < n - K + 1; i++)
        {
            points[i] = KDPoint<T, K>(begin_i + i, begin_i + i + K, 1);
        }
        return points;
    }
}

template<typename T>
T ComputeSum(const vector<T> &data)
{
    unsigned n0 = 1 << 12;
    unsigned p = data.size() / n0;
    T sum = std::accumulate(data.cbegin() + p * n0, data.cend(), 0);
    if (p == 0) return sum;
    else
    {
        vector<T> temp_sum(p, 0);
        for (unsigned i = 0; i < p; i++)
        {
            temp_sum[i] = std::accumulate(
                data.cbegin() + i * n0, data.cbegin() + (i + 1) * n0, 0);
        }
        temp_sum.push_back(sum);
        return ComputeSum(temp_sum);
    }
}

template<typename T>
double ComputeVarience(const vector<T> &data)
{
    vector<double> data_(data.cbegin(), data.cend());
    double avg = ComputeSum(data_) / data.size();
    std::for_each(data_.begin(), data_.end(), [avg](double &x)
    {
        x -= avg; x *= x;
    });
    double sum = ComputeSum(data_);
    sum /= data.size();
    return sum;
}

template<typename T, unsigned K>
void CloseAuxiliaryPoints(vector<KDPoint<T, K> > &points,
                          const vector<unsigned> &rank2index)
{
    const unsigned n = points.size();

    for (unsigned i = 0; i < n; i++)
    {
        if (rank2index[i] >= n - K + 1) points[i].set_count(0);
        else points[i].set_count(1);
    }
}

template<typename T, unsigned K>
void MergeRepeatedPoints(vector<KDPoint<T, K> > &points,
                         const vector<unsigned> &rank2index)
{
    const unsigned n = points.size();

    unsigned k = 0, i = 0;;
    while (i < n)
    {
        int count = 0;
        int count_auxiliary = 0;
        while (k < n && points[k] == points[i])
        {
            // Auxiliary points. 
            if (rank2index[k] >= n - K + 1) count_auxiliary++;
            points[k].set_count(0);
            k++;
        }
        points[k - 1].set_count(k - i - count_auxiliary);
        i = k;
    }
}

template<typename T, unsigned K>
vector<KDPoint<unsigned, K - 1> > Map2Grid(
    const vector<KDPoint<T, K> > &points, const vector<unsigned> &rank2index)
{
    const unsigned n = points.size();

    assert(n == rank2index.size());
    vector<KDPoint<unsigned, K - 1> > result(n);

    vector<unsigned> index2rank = GetInverseMap(rank2index);
    // Mapping q in the paper. 
    vector<unsigned> rank2next(n + 1);
    for (unsigned i = 0; i < n; i++)
    {
        rank2next[i] = index2rank[(rank2index[i] + 1) % n];
    }
    for (unsigned i = 0; i < n; i++)
    {
        if (points[i].count() == 0) continue;
        result[i].set_count(points[i].count());
        unsigned grid = i;
        for (unsigned j = 0; j < K - 1; j++)
        {
            grid = rank2next[grid];
            result[i][j] = grid;
        }
    }
    return result;
}


template<typename T, unsigned K>
Bounds GetRankBounds(const vector<KDPoint<T, K> > &points, T r)
{
    size_t n = points.size();
    Bounds bounds(n);
    vector<T> data(n);
    for (size_t i = 0; i < n; i++)
        data[i] = points[i][0];

    size_t k = 0;
    for (size_t i = 0; i < n; i++)
    {
        while (data[k] + r < data[i]) k++;
        bounds.lower_bounds[i] = k;
    }
    k = n - 1;
    for (size_t i = n; i > 0; i--)
    {
        while (data[k] - r > data[i - 1]) k--;
        bounds.upper_bounds[i - 1] = k;
    }
    return bounds;
}

template<unsigned K>
Range<unsigned, K> GetHyperCube(const KDPoint<unsigned, K> &point,
                                const Bounds &bounds)
{
    Range<unsigned, K> result;
    for (size_t i = 0; i < K; ++i)
    {
        result.lower_ranges[i] = bounds.lower_bounds[point[i]];
        result.upper_ranges[i] = bounds.upper_bounds[point[i]];
    }
    return result;
}


} // namespace kdtree_mddc

#endif // !__UTILS__
