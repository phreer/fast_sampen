# distutils language = c++
# distutils include_dirs = ../include

from libcpp.vector cimport vector
from libcpp cimport bool
from sample_entropy_calculator cimport SampleEntropyCalculatorMao
from sample_entropy_calculator cimport SampleEntropyCalculatorRKD
from sample_entropy_calculator cimport SampleEntropyCalculatorFastDirect
from sample_entropy_calculator cimport SampleEntropyCalculatorDirect
from sample_entropy_calculator cimport SampleEntropyCalculatorSamplingDirect

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

cdef class SampEnSKD:
  cdef SampleEntropyCalculatorMao[double]* c_

  """Compute sample entropy with sliding kd tree method."""
  def __cinit__(self, const vector[double] &data, double r, unsigned m,
                OutputLevel level):
    self.c_ = new SampleEntropyCalculatorMao[double](data, r, m, level)

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

cdef class SampEnRKD:
  cdef SampleEntropyCalculatorRKD[double]* c_

  """Compute sample entropy with range-kd tree method."""
  def __cinit__(self, const vector[double] &data, double r, unsigned m,
                OutputLevel level):
    self.c_ = new SampleEntropyCalculatorRKD[double](data, r, m, level)
  
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

cdef class SampEnFD:
  cdef SampleEntropyCalculatorFastDirect[double]* c_

  """Compute sample entropy with fast direct method."""
  def __cinit__(self, const vector[double] &data, double r, unsigned m,
                OutputLevel level):
    self.c_ = new SampleEntropyCalculatorFastDirect[double](data, r, m, level)
  
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

cdef class SampEnD:
  cdef SampleEntropyCalculatorDirect[double]* c_

  """Compute sample entropy with (trivial) direct method."""
  def __cinit__(self, const vector[double] &data, double r, unsigned m,
                OutputLevel level):
    self.c_ = new SampleEntropyCalculatorDirect[double](data, r, m, level)
  
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


cdef class SampEnSamplingD:
  cdef SampleEntropyCalculatorSamplingDirect[double]* c_

  """Compute sample entropy with sampling and (trivial) direct method."""
  def __cinit__(self, const vector[double] &data, double r, unsigned m,
                unsigned sample_size, unsigned sample_num,
                double real_entropy, double real_a_norm, double real_b_norm,
                RandomType random_type, bool random_, bool presort,
                OutputLevel level):
    self.c_ = new SampleEntropyCalculatorSamplingDirect[double](
        data, r, m, sample_size, sample_num,
        real_entropy, real_a_norm, real_b_norm,
        random_type, random_, presort, level)
  
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
  
  def a_list(self):
    return self.c_.get_a_vec()

  def b(self):
    return self.c_.get_b()
    
  def b_norm(self):
    return self.c_.get_b_norm()

  def b_list(self):
    return self.c_.get_b_vec()