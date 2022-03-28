# distutils language = c++
# distutils include_dirs = ../include

from libcpp.vector cimport vector
from libcpp cimport bool
from sample_entropy_calculator2d cimport OutputLevel
from sample_entropy_calculator2d cimport SampleEntropyCalculator2DDirect
from sample_entropy_calculator2d cimport SampleEntropyCalculator2DSamplingDirect


cdef extern from "utils.h" namespace "sampen":
  cpdef enum OutputLevel 'sampen::OutputLevel':
    Silent,
    Info,
    Debug


cdef class SampEn2DD:
  cdef SampleEntropyCalculator2DDirect[int]* c_

  """Compute sample entropy (2D) with direct method."""
  def __cinit__(self, const vector[int] &data, int r, unsigned m,
                unsigned width, unsigned height, unsigned moving_step_size,
                unsigned dilation_factor, OutputLevel level):
    self.c_ = new SampleEntropyCalculator2DDirect[int](
        data, r, m, width, height, moving_step_size, dilation_factor, level)
  
  def __dealloc__(self):
    if self.c_ != NULL:
      del self.c_  

  def method_name(self):
    return self.c_.get_method_name()

  def compute(self):
    self.c_.ComputeSampleEntropy()

  def time(self):
    return self.c_.get_computation_time()
  
  def entropy(self):
    return self.c_.get_entropy()

  def a(self):
    return self.c_.get_a()

  def a_norm(self):
    return self.c_.get_a_norm()

  def b(self):
    return self.c_.get_b()

  def b_norm(self):
    return self.c_.get_b_norm()

cdef class SampEn2DSamplingD:
  cdef SampleEntropyCalculator2DSamplingDirect[int]* c_

  """Estimate sample entropy (2D) with sampling and direct method."""
  def __cinit__(self, const vector[int] &data, int r, unsigned m,
                unsigned width, unsigned height, unsigned moving_step_size,
                unsigned dilation_factor,
                unsigned sample_size, unsigned sample_num,
                double real_entropy, double real_a_norm, double real_b_norm, 
                bool random_, OutputLevel level):
    self.c_ = new SampleEntropyCalculator2DSamplingDirect[int](
        data, r, m, width, height, moving_step_size, dilation_factor,
        sample_size, sample_num, real_entropy, real_a_norm, real_b_norm,
        random_, level)

  def __dealloc__(self):
    if self.c_ != NULL:
      del self.c_  
  
  def method_name(self):
    return self.c_.get_method_name()

  def compute(self):
    self.c_.ComputeSampleEntropy()

  def time(self):
    return self.c_.get_computation_time()
  
  def entropy(self):
    return self.c_.get_entropy()

  def a(self):
    return self.c_.get_a()

  def a_norm(self):
    return self.c_.get_a_norm()

  def b(self):
    return self.c_.get_b()

  def b_norm(self):
    return self.c_.get_b_norm()