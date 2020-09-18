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
#include <sstream>
#include <numeric>
#include <algorithm>
#include <limits>
#include <assert.h>

#include "kdpoint.h"

#ifdef ENABLE_DEBUG_MACRO
#define DEBUG
#include <iostream> 
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
vector<T> ReadData(std::string filename, std::string input_type = "simple", unsigned n = 0);

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
vector<KDPoint<T, K> > GetKDPoints(const vector<T> &data, int count = 1);

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
    const vector<KDPoint<T, K> > &points, const vector<unsigned> &rank2index,
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
template<typename T, unsigned K>
Bounds GetRankBounds(const vector<KDPoint<T, K> > &points, T r);

/*
* @brief Given a point (in grid), get the bound.
*/
template<unsigned K>
Range<unsigned, K> GetHyperCube(
        const KDPoint<unsigned, K> &point, const Bounds &bounds);


class ArgumentParser
{
public:
    ArgumentParser(int argc, char *argv[]) : arg_list(argv, argv + argc) {}
    string getArg(const string &arg);
    bool isOption(const string &opt);
    long getArgLong(const string &arg, long default_); 
    double getArgDouble(const string &arg, double default_); 
private:
    vector<string> arg_list;
};

//////////////////////////////////////////////////////////////////////////
// Implementation 
//////////////////////////////////////////////////////////////////////////

template<typename T>
vector<T> ReadData(std::string filename, std::string input_type, unsigned n, 
                   unsigned line_offset)
{
    ifstream ifs(filename);
    vector<T> result;
    if (!ifs.is_open())
    {
        std::cerr << "Cannot open file! (filename: " << filename << ")\n";
        exit(-1);
    } 
    if (n == 0) n = std::numeric_limits<unsigned>::max();
    unsigned count = 0;
    T x = 0;
    if (input_type == "simple")
    {
        while (count < n + line_offset && ifs >> x) 
        {
            if (count >= line_offset) result.push_back(x);
            ++count;
        }
    } else if (input_type == "multirecord") 
    {
        std::string line;
        while(count < n + line_offset && std::getline(ifs, line))
        {
            std::istringstream iss(line);
            if (!(iss >> x >> x)) 
            {
                std::cerr << "Input file foramt error. \n";
                exit(-1); 
            }
            if (count >= line_offset) result.push_back(x);
            ++count;
        }
    } else 
    {
        cerr << "Invalid argument n. \n"; 
        cerr << "File: " << __FILE__ << ", Line:: " << __LINE__ << std::endl;
    }
    ifs.close();
    return result;
}

template<typename T, unsigned K>
vector<KDPoint<T, K> >GetKDPoints(const vector<T> &data, int count)
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
            points[i] = KDPoint<T, K>(begin_i + i, begin_i + i + K, count);
        }
        return points;
    }
}

template<typename T>
T ComputeSum(const vector<T> &data)
{
    unsigned n0 = 1 << 4;
    unsigned p = data.size() / n0;
    T sum = std::accumulate(data.cbegin() + p * n0, data.cend(), (T) 0);
    if (p == 0) return sum;
    else
    {
        vector<T> temp_sum(p, 0);
        for (unsigned i = 0; i < p; i++)
        {
            temp_sum[i] = std::accumulate(
                data.cbegin() + i * n0, data.cbegin() + (i + 1) * n0, (T) 0);
        }
        temp_sum.push_back(sum);
        return ComputeSum(temp_sum);
    }
}

template<typename T>
double ComputeVarience(const vector<T> &data)
{
    vector<long double> data_(data.cbegin(), data.cend());
    long double avg = ComputeSum(data_) / data.size();
    std::for_each(data_.begin(), data_.end(), [avg](long double &x)
    {
        x -= avg; x *= x;
    });
    long double sum = ComputeSum(data_);
    #ifdef DEBUG
    std::cout << "avg: " << avg << std::endl;
    std::cout << "variance * size: " << sum << std::endl;
    #endif
    sum /= data.size();
    return (double) sum;
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
        int count_auxiliary = 0;
        int count = 0;
        while (k < n && points[k] == points[i])
        {
            // Auxiliary points. 
            if (rank2index[k] >= n - K + 1) count_auxiliary++;
            count += points[k].count(); 
            points[k].set_count(0);
            k++;
        }
        points[k - 1].set_count(count - count_auxiliary);
        i = k;
    }
}

template<typename T, unsigned K>
vector<KDPoint<unsigned, K - 1> > Map2Grid(
    const vector<KDPoint<T, K> > &points, const vector<unsigned> &rank2index,
    bool skip_nocount)
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
    if (skip_nocount)
    {
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
    }
    else
    {
        for (unsigned i = 0; i < n; i++)
        {
            result[i].set_count(points[i].count());
            unsigned grid = i;
            for (unsigned j = 0; j < K - 1; j++)
            {
                grid = rank2next[grid];
                result[i][j] = grid;
            }
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
