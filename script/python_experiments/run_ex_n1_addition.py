from multiprocessing.connection import wait
import os
import sys
import numpy as np
import sampen
from scipy.io import loadmat
from concurrent.futures import ProcessPoolExecutor

from exp_utils import estimate_sampen_and_save_statistics

n_workers = 10
r = 0.15
m = 5
tag = '220325_m%d' % m
script_dir = os.path.dirname(sys.argv[0])
input_dir = os.path.join(os.environ['HOME'], 'workspace/entropy/additional_data')
output_dir = os.path.join(script_dir, '../../result/convergence_n1_addition', tag)
n0 = 2000

output_name = '{db_name}_{record_name}_n1{n1}'

def get_additional_record(db_name, record_name=None):
  input_path = os.path.join(input_dir, '%s.mat' % db_name)
  try:
    d = loadmat(input_path)
  except Exception as e:
    print('Failed to load database %s.' % db_name, file=sys.stderr)
    return None
  if not record_name in d.keys():
    print('Record %s not found in database %s.' % (record_name, db_name), file=sys.stderr)
    return None
  return d[record_name].flatten()
    

def do_convergence_n1_additional(db_name, record_name):
  list_n1 = [i * 10 for i in range(1, 26)]
  d = get_additional_record(db_name, record_name)
  if d is None:
    return
  r_scaled = r * np.std(d)
  s = sampen.SampEnRKD(d, r_scaled, m, 0)
  r_e, r_a, r_b = s.entropy(), s.a_norm(), s.b_norm()
  print('db: %s, record: %s, r_e: %.4f, r_a: %.4f, r_b: %.4f' % (db_name, record_name, r_e, r_a, r_b))
  pool = ProcessPoolExecutor(max_workers=n_workers)
  
  for n1 in list_n1:
    def callback(future):
      if future.exception():
        print('error occurred %s %s %d' % (db_name, record_name, n1))
        print(future)
      else:
        print('completed %s %s %d' % (db_name, record_name, n1))
    output_filename = output_name.format(db_name=db_name, record_name=record_name, n1=n1)
    print('submit %s %s %d' % (db_name, record_name, n1))
    ft = pool.submit(estimate_sampen_and_save_statistics, sig=d, r=r, m=m,
                     sig_offset=0, n0=n0, n1=n1, outputdir=output_dir,
                     output_filename=output_filename, groundtruth=(r_e, r_a, r_b))
    ft.add_done_callback(callback)
  return pool

  
if __name__ == '__main__':
  os.makedirs(output_dir, exist_ok=True)
  list_db_records = [('kmni_climate', 'uurgeg_t'),
                     ('gears', 'Miss_30_2_1'),
                     ('rolling_bearing', 'X110_DE_time')]
  list_pool = []
  for db_name, record_name in list_db_records:
    pool = do_convergence_n1_additional(db_name=db_name, record_name=record_name)
    list_pool.append(pool)