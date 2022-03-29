/* file: experiment_n0n1.cpp
 * date: 2021-01-06
 * author: phree
 *
 * description: This program do experiement on convergence of sample size N_0
 * and sample num N_1.
 */
#include <iostream>
#include <string>
#include <time.h>

#include "experiment.h"
#include "random_sampler.h"
#include "sample_entropy_calculator_direct.h"
#include "sample_entropy_calculator_kd.h"
#include "utils.h"

using namespace sampen;
using std::cout;
using std::endl;
using std::string;
using std::vector;

// clang-format: off
char usage[] =
    "Usage: %s --input <INPUT> --input-type {simple, multirecord}\\\n"
    "                   -r <THRESHOLD> -m <TEMPLATE_LENGTH>\\\n"
    "                   [-output-level {1,2,3}]\\\n"
    "Arguments:\n"
    "--input <INPUT>         The file name of the input file.\n"
    "-n <N>                  The length of signal to read.\n"
    "--input-format <FORMAT> The format of the input file. Should be either simple\n"
    "                        or multirecord. If set to simple, then each line of the\n"
    "                        input file contains exactly one column; if set to\n"
    "                        multirecord, then each line contains <NUM_RECORD> + 1\n"
    "                        columns, of which the first indicates the line number\n"
    "                        and the remaining columns are instances of the records.\n"
    "                        The default value is simple.\n"
    "--input-type <TYPE>     The data type of the input data, either int or float.\n"
    "                        Default: double.\n"
    "-r <R>                  The threshold argument in sample entropy.\n"
    "-m <M>                  The template length argument of sample entropy. Note\n"
    "                        that this program only supports 2 <= m <= 10.\n"
    "--output-level <LEVEL>  The amount of information printed. Should be one of\n"
    "                        {0,1,2}. Level 0 is most silent while level 2 is for\n"
    "                        debugging.\n"
    "--quasi-type <TYPE>     The type of the quasi-random sequence for sampling,\n"
    "                        can be one of the following: sobol, halton,\n"
    "                        reversehalton or niederreiter_2. Default: sobol.\n\n"
    "--computation-times COMPUTATION_TIMES\n"
    "    The number of computations for mean and variance of the sampling methods.\n"
    "--sample-size-array N01,N02,N03,N0n\n"
    "    An array of sample sizes (N_0s) that computations should take on. "
    "The\n"
    "    default value is `200,400,600,...,4000`. Note that your command should not\n"
    "    contain `...`.\n"
    "--sample-num-array N11,N12,N13,N1n\n"
    "    An array of sample numbers (N_1) that computations should take on. The\n"
    "    default value is `10,20,...,250`. Note that your command should not contain\n"
    "    `...`.\n"
    "Options:\n"
    "--random\n"
    "    If this option is enabled, the random seed will be set randomly.\n"
    "-q\n"
    "    If this option is enabled, the quasi-Monte Carlo based method will be\n"
    "    conducted.\n"
    "-u | --uniform\n"
    "    If this option is enabled, the Monte Carlo based method using uniform\n"
    "    distribution will be conducted.\n"
    "--kdtree-sample\n"
    "    If this option is enabled, the kd tree based sampling method will be\n"
    "    conducted.\n"
    "--swr\n"
    "    If this option is enabled, then the Monte Carlo based method without\n"
    "    placement will be performed.\n"
    "--grid\n"
    "    If this option is enabled, then the quasi-Monte Carlo based method using\n"
    "    grid (lattice) as sampling indexes will be performed.\n"
    "--presort\n"
    "    If this option is enabled, then a presorting operation is conducted\n"
    "    before sampling in quasi-Monte Carlo based method.\n"
    "--variance\n"
    "    If this option is enabled, then the variance of the results of sampling\n"
    "    methods will be computed.\n"
    "--help | -h\n"
    "    Show this message.\n";

template <typename T>
void PrintSampenSetting(unsigned line_offset, unsigned n, T r, unsigned K,
                        std::string filename, double var);

void ParseArgument(int argc, char *argv[]);

struct Argument {
  unsigned template_length;
  unsigned line_offset;
  string filename;
  string input_format;
  string input_type;
  unsigned data_length;
  vector<unsigned> sample_sizes;
  vector<unsigned> sample_nums;
  double r;
  OutputLevel output_level;
  bool random_, variance;
  bool q, u, swr, presort, grid;
  bool kdtree_sample;
  RandomType rtype;
  void PrintArguments() const;
} arg;

template <typename T> void SampleEntropyN0N1();

int main(int argc, char *argv[]) {
#ifdef DEBUG
  cout << "This is a debug version. " << std::endl;
#endif
  ParseArgument(argc, argv);

  if (arg.input_type == "double") {
    using Type = double;
    SampleEntropyN0N1<Type>();
  } else if (arg.input_type == "int") {
    using Type = int;
    SampleEntropyN0N1<Type>();
  } else {
    cerr << "Invalid argument: -type " << arg.input_type << ".\n";
    cerr << "File: " << __FILE__ << ", Line: " << __LINE__ << std::endl;
    exit(-1);
  }
  return 0;
}

void Argument::PrintArguments() const {
  std::cout << "Sample Entropy Computation Setting: " << std::endl;
  std::cout << "\tfilename: " << arg.filename << std::endl;
  std::cout << "\tinput type: " << arg.input_type << std::endl;
  std::cout << "\tline offset: " << arg.line_offset << std::endl;
  std::cout << "\tdata length: " << arg.data_length << std::endl;
  std::cout << "\ttemplate length: " << arg.template_length << std::endl;
  std::cout << "\tthreshold: " << arg.r << std::endl;
  std::cout << "\trandom: " << arg.random_ << std::endl;
  std::cout << "\tquasi type: " << random_type_names[arg.rtype] << std::endl;
  if (!arg.sample_sizes.empty()) {
    std::cout << "\tsample size array: ";
    for (auto sample_size : sample_sizes) {
      std::cout << sample_size << " ";
    }
    std::cout << "\n";
  }
  if (!arg.sample_nums.empty()) {
    std::cout << "\tsample num array: ";
    for (auto sample_nums : sample_nums) {
      std::cout << sample_nums << " ";
    }
    std::cout << "\n";
  }
}

void ParseArgument(int argc, char *argv[]) {
  ArgumentParser parser(argc, argv);
  char _usage[sizeof(usage) + 256];
  sprintf(_usage, usage, argv[0]);
  if (parser.isOption("--help") || parser.isOption("-h")) {
    std::cout << _usage;
    exit(0);
  }

  arg.filename = parser.getArg("--input");
  if (arg.filename.size() == 0) {
    cerr << "Specify a filename with --input <INPUT>." << endl;
    cerr << _usage;
    exit(-1);
  }

  arg.input_format = parser.getArg("--input-format");
  if (arg.input_format.size() == 0)
    arg.input_format = "simple";
  if (arg.input_format != "simple" && arg.input_format != "multirecord") {
    cerr << "Invalid argument: --input-format " << arg.input_format;
    cerr << ", should be either simple or multirecord. \n";
    exit(-1);
  }

  arg.input_type = parser.getArg("--input-type");
  if (arg.input_type.size() == 0)
    arg.input_type = "double";

  arg.data_length = parser.getArgLong("-n", 1000010);

  arg.r = parser.getArgDouble("-r", -1.);
  if (arg.r < 0) {
    cerr << "Specify a positive threshold with -r <R>. " << endl;
    cerr << _usage;
    exit(-1);
  }

  string template_length = parser.getArg("-m");
  if (template_length.size() == 0) {
    cerr << "Specify template length with -m <M> (1 <= M <= 10." << endl;
    cerr << _usage;
    exit(-1);
  } else {
    arg.template_length = static_cast<unsigned>(std::stoi(template_length));
    if (arg.template_length < 1 || arg.template_length > 10) {
      cerr << "Specify template length with -m <M> (1 <= M <= 10." << endl;
      exit(-1);
    }
  }

  arg.line_offset =
      static_cast<unsigned>(parser.getArgLong("--line-offset", 0));

  string output_level = parser.getArg("--output-level");
  if (output_level.size() == 0)
    arg.output_level = Info;
  else if (output_level == "0")
    arg.output_level = Silent;
  else if (output_level == "1")
    arg.output_level = Info;
  else if (output_level == "2")
    arg.output_level = Debug;
  else {
    cerr << "Invalid argument --output-level " << output_level;
    cerr << ", must be one of the following: {1, 2, 3}. \n";
    cerr << _usage;
    exit(-1);
  }

  arg.q = parser.isOption("-q");
  arg.u = parser.isOption("-u") || parser.isOption("--uniform");
  arg.swr = parser.isOption("--swr");
  arg.kdtree_sample = parser.isOption("--kdtree-sample");
  arg.grid = parser.isOption("--grid");
  if (arg.q || arg.u || arg.swr || arg.grid || arg.kdtree_sample) {
    vector<int> sample_sizes = parser.getArgIntArray("--sample-size-array");
    if (sample_sizes.empty()) {
      for (unsigned i = 1; i <= 20; ++i) {
        arg.sample_sizes.push_back(i * 200);
      }
    } else {
      for (auto x : sample_sizes) {
        if (x <= 0) {
          std::cerr << "Invalid sample size (N0): " << x << std::endl;
          exit(-1);
        }
        arg.sample_sizes.push_back(static_cast<unsigned>(x));
      }
    }
    vector<int> sample_nums = parser.getArgIntArray("--sample-num-array");
    if (sample_nums.empty()) {
      for (unsigned i = 1; i <= 25; ++i) {
        arg.sample_nums.push_back(i * 10);
      }
    } else {
      for (auto x : sample_nums) {
        if (x <= 0) {
          std::cerr << "Invalid sample num (N1): " << x << std::endl;
          exit(-1);
        }
        arg.sample_nums.push_back(static_cast<unsigned>(x));
      }
    }

    arg.random_ = parser.isOption("--random");
    arg.variance = parser.isOption("--variance");
    if (arg.q || arg.grid) {
      arg.presort = parser.isOption("--presort");
      std::string rtype = parser.getArg("--quasi-type");
      if (rtype.size() == 0 || rtype == "sobol")
        arg.rtype = SOBOL;
      else if (rtype == "halton")
        arg.rtype = HALTON;
      else if (rtype == "reverse_halton")
        arg.rtype = REVERSE_HALTON;
      else if (rtype == "niederreiter_2")
        arg.rtype = NIEDERREITER_2;
      else {
        cerr << "Invalid argument --quasi-random " << rtype << ". ";
        cerr << "Should be one of the following: sobol, halton, "
                "reverse_halton, niederreiter_2 or grid. \n";
        exit(-1);
      }
    }
  }
}

template <typename T> void SampleEntropyN0N1() {
  const unsigned K = arg.template_length;
  vector<T> data;
  ReadData<T>(data, arg.filename, arg.input_format, arg.data_length,
              arg.line_offset);

  unsigned n = data.size();
  if (n <= K) {
    std::cerr << "Data length n " << n << " is too short (K = " << K;
    std::cerr << "). \n";
    exit(-1);
  }

  double var = ComputeVariance(data);
  T r_scaled = static_cast<T>(sqrt(var) * arg.r);

  if (r_scaled < 0) {
    cerr << "Invalid r [scaled]: " << r_scaled << ".\n";
    exit(-1);
  }

  cout.precision(4);
  cout << std::scientific;
  cout << "========================================";
  cout << "========================================\n";
  arg.PrintArguments();
  std::cout << "\tvariance: " << var << std::endl;
  std::cout << "\tr (scaled): " << r_scaled << std::endl;

  // Compute sample entropy.
  SampleEntropyCalculatorMao<T> sec(data, r_scaled, K, arg.output_level);
  sec.ComputeSampleEntropy();
  cout << sec.get_result_str();

  unsigned n_computation = 1;
  if (arg.variance)
    n_computation = 50;
  for (unsigned i = 0; i < arg.sample_sizes.size(); ++i) {
    unsigned sample_size = arg.sample_sizes[i];
    for (unsigned j = 0; j < arg.sample_nums.size(); ++j) {
      unsigned sample_num = arg.sample_nums[j];
      cout << "========================================";
      cout << "========================================\n";
      cout << "sample_size: " << sample_size << endl;
      cout << "sample_num: " << sample_num << endl;
      cout << "========================================";
      cout << "========================================\n";
      if (arg.kdtree_sample) {
        SampleEntropyCalculatorSamplingMao<T> secds(
            data, r_scaled, K, sample_size, sample_num,
            sec.get_entropy(), sec.get_a_norm(), sec.get_b_norm(), SWR_UNIFORM,
            arg.random_, arg.output_level);
        SampleEntropySamplingExperiment(secds, n_computation);
      }
      if (arg.u) {
        SampleEntropyCalculatorSamplingDirect<T> secds(
            data, r_scaled, K, sample_size, sample_num,
            sec.get_entropy(), sec.get_a_norm(), sec.get_b_norm(), UNIFORM,
            arg.random_, false, arg.output_level);
        SampleEntropySamplingExperiment(secds, n_computation);
      }

      if (arg.swr) {
        SampleEntropyCalculatorSamplingDirect<T> secds(
            data, r_scaled, K, sample_size, sample_num,
            sec.get_entropy(), sec.get_a_norm(), sec.get_b_norm(), SWR_UNIFORM,
            arg.random_, false, arg.output_level);
        SampleEntropySamplingExperiment(secds, n_computation);
      }
      if (arg.q) {
        SampleEntropyCalculatorSamplingDirect<T> secds(
            data, r_scaled, K, sample_size, sample_num,
            sec.get_entropy(), sec.get_a_norm(), sec.get_b_norm(), arg.rtype,
            arg.random_, false, arg.output_level);

        SampleEntropySamplingExperiment(secds, n_computation);
        if (arg.presort) {
          SampleEntropyCalculatorSamplingDirect<T> secds(
              data, r_scaled, K, sample_size, sample_num,
              sec.get_entropy(), sec.get_a_norm(), sec.get_b_norm(), arg.rtype,
              arg.random_, true, arg.output_level);
          SampleEntropySamplingExperiment(secds, n_computation);
        }
      }

      if (arg.grid) {
        SampleEntropyCalculatorSamplingDirect<T> secds(
            data, r_scaled, K, sample_size, sample_num,
            sec.get_entropy(), sec.get_a_norm(), sec.get_b_norm(), GRID,
            arg.random_, false, arg.output_level);
        SampleEntropySamplingExperiment(secds, n_computation);
        if (arg.presort) {
          SampleEntropyCalculatorSamplingDirect<T> secds(
              data, r_scaled, K, sample_size, sample_num,
              sec.get_entropy(), sec.get_a_norm(), sec.get_b_norm(), GRID,
              arg.random_, true, arg.output_level);
          SampleEntropySamplingExperiment(secds, n_computation);
        }
      }
    }
  }
  cout << "========================================";
  cout << "========================================\n";
}
