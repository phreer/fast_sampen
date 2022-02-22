/**
 * @file utils.cpp
 * @author Weifeng Liu
 * @date 2020/05/09
 *
 * @brief Some utility function and class.
 *
 * @details
 */
#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdexcept>

#include "utils.h"

namespace sampen {
using namespace std;

double ComputeSampen(double A, double B, unsigned N, unsigned m) {
  if (A > 0 && B > 0) {
    return -log(A / B);
  } else {
    return -log((N - m - 1) / (N - m));
  }
}

vector<unsigned> GetInverseMap(const vector<unsigned> &map) {
  vector<unsigned> inv(map);
  size_t n = map.size();
  for (size_t i = 0; i < n; i++) {
    inv[map[i]] = i;
  }
  return inv;
}

string ArgumentParser::getArg(const string &arg) {
  auto iter = std::find(arg_list.cbegin(), arg_list.cend(), arg);
  if (iter != arg_list.cend() && ++iter != arg_list.cend()) {
    return *iter;
  } else
    return string("");
}

bool ArgumentParser::isOption(const string &opt) {
  return find(arg_list.cbegin(), arg_list.cend(), opt) != arg_list.cend();
}

long ArgumentParser::getArgLong(const string &arg, long default_) {
  long result = default_;
  std::string arg_v = getArg(arg);
  if (arg_v.size() == 0)
    return result;
  try {
    result = std::stol(arg_v);
  } catch (const std::invalid_argument &e) {
    std::cerr << "Invalid argument: " << arg << " " << arg_v;
    std::cerr << "\n" << e.what() << '\n';
    exit(-1);
  }
  return result;
}

double ArgumentParser::getArgDouble(const string &arg, double default_) {
  double result = default_;
  std::string arg_v = getArg(arg);
  if (arg_v.size() == 0)
    return result;
  try {
    result = std::stod(arg_v);
  } catch (const std::invalid_argument &e) {
    std::cerr << "Invalid argument: " << arg << " " << arg_v;
    std::cerr << "\n" << e.what() << "\n";
    exit(-1);
  }
  return result;
}

vector<int> ParseIntArrayFromString(const string &arg) {
  vector<int> result;
  if (arg.empty())
    return result;
  int curr = 0;
  int next = std::find(arg.cbegin(), arg.cend(), ',') - arg.cbegin();
  while (true) {
    int val = std::stoi(arg.substr(curr, next - curr));
    result.push_back(val);
    curr = next + 1;
    if (curr == static_cast<int>(arg.size() + 1))
      break;
    next = std::find(arg.cbegin() + next + 1, arg.cend(), ',') - arg.cbegin();
  }
  return result;
}

vector<int> ArgumentParser::getArgIntArray(const string &arg) {
  std::string arg_v = getArg(arg);
  vector<int> result;
  try {
    result = ParseIntArrayFromString(arg_v);
  } catch (const std::invalid_argument &e) {
    std::cerr << "Invalid argument: " << arg << " " << arg_v;
    std::cerr << "\n" << e.what() << "\n";
    exit(-1);
  }
  return result;
}

void PrintSeperator(char x) {
  const int kCount = 80;
  for (unsigned i = 0; i < kCount; ++i) {
    std::cout << x;
  }
  std::cout << std::endl;
}
} // namespace sampen
