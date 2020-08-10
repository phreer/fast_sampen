/**
 * @file utils.cpp
 * @author Weifeng Liu
 * @date 2020/05/09
 * 
 * @brief Some utility function and class. 
 *
 * @details 
 */
#include <iostream>
#include <fstream>
#include <algorithm>

#include "utils.h"

namespace kdtree_mddc
{
using namespace std;

double ComputeSampen(double A, double B, unsigned N, unsigned m, OutputLevel output_level)
{
    if (output_level) std::cout << "A: " << A << ", B: " << B << std::endl;
    if (A > 0 && B > 0)
    {
        return -log(A / B);
    } else
    {
        return -log((N - m - 1) / (N - m));
    }
}

vector<unsigned> GetInverseMap(const vector<unsigned> &map)
{
    vector<unsigned> inv(map);
    size_t n = map.size();
    for (size_t i = 0; i < n; i++)
    {
        inv[map[i]] = i;
    }
    return inv;
}

string ArgumentParser::getArg(const string &arg)
{
    auto iter = std::find(arg_list.cbegin(), arg_list.cend(), arg);
    if (iter != arg_list.cend() && ++iter != arg_list.cend())
    {
        return *iter;
    } else return string("");
}

bool ArgumentParser::isOption(const string &opt)
{
    return find(arg_list.cbegin(), arg_list.cend(), opt) != arg_list.cend();
}

} // namespace kdtree_mddc