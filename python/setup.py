from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import platform

description = 'Fast computation (or estimation) of sample entropy.'
extra_compile_args = ['-std=c++11']
if platform.system() == 'Darwin':
    extra_compile_args += ['-mmacosx-version-min=10.7', '-stdlib=libc++']

setup(name='sampen',
      description=description,
      version='1.0.0',
      author='Phree Liu',
      author_email='iphreeliu@gmail.com',
      ext_modules=cythonize([Extension('sampen',
                                       sources=[
                                         'sample_entropy_calculator.pyx',
                                       ],
                                       language='c++',
                                       include_dirs=['../include'],
                                       library_dirs=['../build/lib', '/usr/lib'],
                                       libraries=['sampen', 'gsl', 'gslcblas'],
                                       extra_compile_args=extra_compile_args),
                             Extension('sampen2d',
                                       sources=[
                                         'sample_entropy_calculator2d.pyx',
                                       ],
                                       language='c++',
                                       include_dirs=['../include'],
                                       library_dirs=['../build/lib', '/usr/lib'],
                                       libraries=['sampen', 'gsl', 'gslcblas'],
                                       extra_compile_args=extra_compile_args),
]))