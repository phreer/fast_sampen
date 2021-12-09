#ifndef __KDTREE_MDDC_EXPERIMENT_H__
#define __KDTREE_MDDC_EXPERIMENT_H__
#include <math.h>

#include "sample_entropy_calculator.h"

namespace kdtree_mddc {

template <typename T, unsigned K>
void SampleEntropySamplingExperiment(
    SampleEntropyCalculatorSampling<T, K> &secds, unsigned n_computation) {
  vector<double> errs_sampen(n_computation);
  vector<double> errs_a(n_computation);
  vector<double> errs_b(n_computation);
  vector<double> computation_times(n_computation);
  for (unsigned i = 0; i < n_computation; ++i) {
    secds.ComputeSampleEntropy();
    if (n_computation == 1)
      std::cout << secds.get_result_str();
    errs_sampen[i] = secds.get_err_entropy();
    errs_a[i] = secds.get_err_a();
    errs_b[i] = secds.get_err_b();
    computation_times[i] = secds.get_computation_time();
  }

  if (n_computation > 1) {
    std::cout << "----------------------------------------"
              << "----------------------------------------\n"
              << secds.get_method_name() << std::endl;
    std::cout << "\terrs_entropy: ";
    for (unsigned i = 0; i < n_computation; ++i) {
      std::cout << errs_sampen[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "\terrs_a: ";
    for (unsigned i = 0; i < n_computation; ++i) {
      std::cout << errs_a[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "\terrs_b: ";
    for (unsigned i = 0; i < n_computation; ++i) {
      std::cout << errs_b[i] << ", ";
    }
    std::cout << std::endl;

    std::cout << "\tcomputation_times: ";
    for (unsigned i = 0; i < n_computation; ++i) {
      std::cout << computation_times[i] << ", ";
    }
    std::cout << std::endl;
    for (unsigned i = 0; i < n_computation; ++i) {
      errs_sampen[i] = fabs(errs_sampen[i]);
      errs_a[i] = fabs(errs_a[i]);
      errs_b[i] = fabs(errs_b[i]);
    }
    double var_errs_sampen = ComputeVariance<double>(errs_sampen);
    double mean_errs_sampen = ComputeSum<double>(errs_sampen) / n_computation;
    double var_errs_a = ComputeVariance<double>(errs_a);
    double mean_errs_a = ComputeSum<double>(errs_a) / n_computation;
    double var_errs_b = ComputeVariance<double>(errs_b);
    double mean_errs_b = ComputeSum<double>(errs_b) / n_computation;
    double mean_computation_time =
        ComputeSum<double>(computation_times) / n_computation;
    std::cout << "\tmean_errs_sampen: " << mean_errs_sampen << std::endl;
    std::cout << "\tstd_errs_sampen: " << sqrt(var_errs_sampen) << std::endl;
    std::cout << "\tmean_errs_a: " << mean_errs_a << std::endl;
    std::cout << "\tstd_errs_a: " << sqrt(var_errs_a) << std::endl;
    std::cout << "\tmean_errs_b: " << mean_errs_b << std::endl;
    std::cout << "\tstd_errs_b: " << sqrt(var_errs_b) << std::endl;
    std::cout << "\tmean_computation_time: " << mean_computation_time
              << std::endl;
    std::cout << "----------------------------------------"
              << "----------------------------------------\n";
  }
}

}; // namespace kdtree_mddc

#endif // __KDTREE_MDDC_EXPERIMENT_H_