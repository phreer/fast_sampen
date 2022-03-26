import os
import numpy as np
import sys
from scipy.io import loadmat
from concurrent.futures import ProcessPoolExecutor

from exp_utils import estimate_sampen_and_save_statistics

if __name__ == '__main__':
  r = 0.15
  n_workers = 30
  list_m = [2, 3, 4, 5]
  n0 = 2048
  n1 = 150
  script_dir = os.path.dirname(sys.argv[0])
  input_dir = os.path.join(os.environ['HOME'], 'workspace/entropy/')
  output_dir = os.path.join(script_dir, '../../result/mix_p_err/n0%d_n1%d' % (n0, n1))
  os.makedirs(output_dir, exist_ok=True)
  pool = ProcessPoolExecutor(max_workers=n_workers)
  list_future = []
  for i in range(21):
    p = 0.05 * i
    d = loadmat(os.path.join(input_dir, 'max_p_%.2f.mat' % p))['mix_p']
    for m in list_m:
      output_filename = 'p%.2f_m%d_r%.2f_n0%d_n1%d.txt' % (p, m, r, n0, n1)
      future = pool.submit(estimate_sampen_and_save_statistics, sig=d, r=r, m=m,
                           sig_offset=0, n0=n0, n1=n1, outputdir=output_dir,
                           output_filename=output_filename)
      list_future.append(future)
      print('submit %s' % output_filename)
  for future in list_future:
    exception = future.exception()
    if exception:
      print(exception)