#include <algorithm>
#include <chrono>
#include <iostream>

#include <cstdlib>
#include <gsl/gsl_qrng.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "random_sampler.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

int main() {
  unsigned pop_size = 1000;
  unsigned sample_size = 20;
  unsigned sample_num = 4;
  auto samples = GetSampleIndicesWR(GRID, pop_size, sample_size, sample_num);
  cout << "pop_size:" << samples.size() << endl;
  cout << "size:" << samples.size() << endl;
  cout << "results: \n";
  for (unsigned i = 0; i < sample_num; i++) {
    cout << "sample " << i << ": ";
    for (unsigned j = 0; j < sample_size; ++j) {
      if (j)
        cout << ", ";
      cout << samples[i][j];
    }
    cout << endl;
  }
}