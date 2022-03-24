# distutils: language = c++
# distutils include_dirs = ../include
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "random_sampler.h":
  cpdef enum RandomType 'RandomType':
    UNIFORM,
    SOBOL,
    HALTON,
    REVERSE_HALTON,
    NIEDERREITER_2,
    GRID,
    SWR_UNIFORM

cdef extern from "utils.h" namespace "sampen":
  cpdef enum OutputLevel 'sampen::OutputLevel':
    Silent,
    Info,
    Debug

cdef extern from "sample_entropy_calculator.h" namespace "sampen":
  cdef cppclass SampleEntropyCalculator[double]:
    pass
  cdef cppclass SampleEntropyCalculatorSampling[double](SampleEntropyCalculator[double]):
    pass

cdef extern from "sample_entropy_calculator_kd.h" namespace "sampen":
  cdef cppclass SampleEntropyCalculatorMao[double](SampleEntropyCalculator[double]):
    SampleEntropyCalculatorMao(const vector[double] &, double, unsigned, OutputLevel) except +
    double get_computation_time()
    double get_entropy()
    long long get_a()
    long long get_b()
    double get_a_norm()
    double get_b_norm()
    string get_method_name()
    void ComputeSampleEntropy()

  cdef cppclass SampleEntropyCalculatorRKD[double](SampleEntropyCalculator[double]):
    SampleEntropyCalculatorRKD(const vector[double] &, double, unsigned, OutputLevel) except +
    double get_computation_time()
    double get_entropy()
    long long get_a()
    long long get_b()
    double get_a_norm()
    double get_b_norm()
    string get_method_name()
    void ComputeSampleEntropy()

cdef extern from "sample_entropy_calculator_direct.h" namespace "sampen":
  cdef cppclass SampleEntropyCalculatorFastDirect[double](SampleEntropyCalculator[double]):
    SampleEntropyCalculatorFastDirect(const vector[double] &, double, unsigned, OutputLevel) except +
    double get_computation_time()
    double get_entropy()
    long long get_a()
    long long get_b()
    double get_a_norm()
    double get_b_norm()
    string get_method_name()
    void ComputeSampleEntropy()

  cdef cppclass SampleEntropyCalculatorDirect[double](SampleEntropyCalculator[double]):
    SampleEntropyCalculatorDirect(const vector[double] &, double, unsigned, OutputLevel) except +
    double get_computation_time()
    double get_entropy()
    long long get_a()
    long long get_b()
    double get_a_norm()
    double get_b_norm()
    string get_method_name()
    void ComputeSampleEntropy()

  cdef cppclass SampleEntropyCalculatorDirect[double](SampleEntropyCalculator[double]):
    SampleEntropyCalculatorDirect(const vector[double] &, double, unsigned, OutputLevel) except +
    double get_computation_time()
    double get_entropy()
    long long get_a()
    long long get_b()
    double get_a_norm()
    double get_b_norm()
    string get_method_name()
    void ComputeSampleEntropy()

  cdef cppclass SampleEntropyCalculatorSamplingDirect[double](SampleEntropyCalculatorSampling[double]):
    SampleEntropyCalculatorSamplingDirect(
      const vector[double] &data, double r, unsigned m,
      unsigned sample_size, unsigned sample_num,
      double real_entropy, double real_a_norm,
      double real_b_norm, RandomType rtype, bool random_, bool presort,
      OutputLevel output_level) except +
    double get_computation_time()
    double get_entropy()
    long long get_a()
    long long get_b()
    double get_a_norm()
    double get_b_norm()
    vector[long long] get_a_vec()
    vector[long long] get_b_vec()
    string get_method_name()
    void ComputeSampleEntropy()
