import os
from os.path import join as pjoin
import numpy as np
import sys
from scipy.io import loadmat
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

from exp_utils import estimate_sampen_and_save_statistics
from exp_utils import read_record

if __name__ == '__main__':
  r = 0.15
  n_workers = 30
  list_m = [2, 3, 4, 5]
  n0 = 2048
  n1 = 150
  n = 1024 * 1024
  script_dir = os.path.dirname(sys.argv[0])
  input_dir = os.path.join(script_dir, '../../data.PhysioNet')
  output_dir = os.path.join(script_dir, '../../result/r_std/n0%d_n1%d' % (n0, n1))
  os.makedirs(output_dir, exist_ok=True)
  pool = ProcessPoolExecutor(max_workers=n_workers)
  list_future = []
  for i in range(21):
    r = 0.05 * i
    input_filename = pjoin(input_dir, 'gaussian/gaussian-1m.txt')
    d = read_record(input_filename, input_format='multirecord', n=n)
    for m in list_m:
      output_filename = 'gaussian_m%d_r%.2f_n0%d_n1%d.txt' % (m, r, n0, n1)
      future = pool.submit(estimate_sampen_and_save_statistics, sig=d, r=r, m=m,
                           sig_offset=0, n0=n0, n1=n1, outputdir=output_dir,
                           output_filename=output_filename)
      list_future.append(future)
      print('submit %s' % output_filename)
  for future in list_future:
    exception = future.exception()
    if exception:
      print(exception)