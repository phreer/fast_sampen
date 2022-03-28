
# distutils: language = c++
# distutils include_dirs = ../include
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "utils.h" namespace "sampen":
  cpdef enum OutputLevel 'sampen::OutputLevel':
    Silent,
    Info,
    Debug

cdef extern from "sample_entropy_calculator2d.h" namespace "sampen":
  cdef cppclass SampleEntropyCalculator2D[T]:
    pass
  cdef cppclass SampleEntropyCalculator2DSampling[T](SampleEntropyCalculator2D[T]):
    pass

  cdef cppclass SampleEntropyCalculator2DDirect[T](SampleEntropyCalculator2D[T]):
    SampleEntropyCalculator2DDirect(const vector[T] &data,
                                    T r, unsigned m,
                                    unsigned width, unsigned height,
                                    unsigned moving_step_size,
                                    unsigned dilation_factor,
                                    OutputLevel output_level) except +
    double get_computation_time()
    double get_entropy()
    long long get_a()
    long long get_b()
    double get_a_norm()
    double get_b_norm()
    string get_method_name()
    void ComputeSampleEntropy()

  cdef cppclass SampleEntropyCalculator2DSamplingDirect[T](SampleEntropyCalculator2DSampling[T]):
    SampleEntropyCalculator2DSamplingDirect(const vector[T] &data,
                                            T r, unsigned m,
                                            unsigned width,
                                            unsigned height,
                                            unsigned moving_step_size,
                                            unsigned dilation_factor,
                                            unsigned sample_size,
                                            unsigned sample_num,
                                            double real_entropy,
                                            double real_a_norm,
                                            double real_b_norm,
                                            bool random_,
                                            OutputLevel output_level) except +
    double get_computation_time()
    double get_entropy()
    long long get_a()
    long long get_b()
    double get_a_norm()
    double get_b_norm()
    string get_method_name()
    void ComputeSampleEntropy()