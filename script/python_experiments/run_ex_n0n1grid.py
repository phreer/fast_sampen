import imp
import os
import sys
from time import sleep
import numpy as np
import sampen
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import Future

from exp_utils import RecordType, estimate_sampen_and_save_statistics
from exp_utils import all_db_records, read_record

n_workers = 36
r = 0.15
script_dir = os.path.dirname(sys.argv[0])
n = 1024 * 1024

def task_groundtruth(db_name, record_name, record_type, m):
  d = read_record(db_name, record_name, record_type, n)
  if d is None:
    raise ValueError('failed to load %s %s' % (db_name, record_name))
  r_scaled = r * np.std(d)
  s = sampen.SampEnRKD(d, r_scaled, m, 0)
  r_e, r_a, r_b = s.entropy(), s.a_norm(), s.b_norm()
  print('db: %s, record: %s, r_e: %.4f, r_a: %.4f, r_b: %.4f' %
      (db_name, record_name, r_e, r_a, r_b))
  return (r_e, r_a, r_b), (db_name, record_name, record_type)

def task_n0n1grid(db_name, record_name, record_type, m, n0, n1, groundtruth):
  d = read_record(db_name, record_name, record_type, n)
  if d is None:
    raise ValueError('failed to load %s %s' % (db_name, record_name))
  output_name = '{db_name}_{record_name}_n0{n0}_n1{n1}'
  output_filename = output_name.format(
      db_name=db_name, record_name=record_name, n0=n0, n1=n1)
  print('submit %s %s m: %d n0: %d n1: %d' % (db_name, record_name, m, n0, n1))
  estimate_sampen_and_save_statistics(
      sig=d, r=r, m=m, sig_offset=0, n0=n0, n1=n1, outputdir=output_dir,
      output_filename=output_filename, groundtruth=groundtruth)
  return (db_name, record_name, m, n0, n1)

def callback_n0n1grid(ft: Future):
  if ft.exception():
    print(ft.exception())
  else:
    db_name, record_name, m, n0, n1 = ft.result()
    print('completed task grid %s %s m %d n0 %d n1 %d' % (db_name, record_name, m, n0, n1))
    

def submit_n0n1grid_tasks(pool, db_name, record_name, record_type, m, groundtruth):
  list_n0 = [i * 200 for i in range(1, 21)]
  list_n1 = [i * 10 for i in range(1, 26)]
  for n1 in list_n1:
    for n0 in list_n0:
      ft = pool.submit(task_n0n1grid, db_name=db_name, record_name=record_name,
                       m=m, record_type=record_type, n0=n0, n1=n1,
                       groundtruth=groundtruth)
      ft.add_done_callback(callback_n0n1grid)
  

def callback_groundtruth(ft: Future):
  if ft.exception():
    print(ft.exception())
  else:
    groundtruth = ft.result()[0]
    db_name, record_name, record_type = ft.result()[1]
    if groundtruth is None:
      print('failed to get groundtruth', file=sys.stderr)
      pool.shutdown(wait=False)
      exit(1)
    submit_n0n1grid_tasks(pool, db_name, record_name, record_type, m, groundtruth)
    print('completed task groundtruth %s %s %d' % (db_name, record_name, m))

if __name__ == '__main__':
  pool = ProcessPoolExecutor(max_workers=n_workers)
  list_m = range(4, 5)
  list_futures = []
  for m in list_m:
    tag = '220328_n%d_m%d' % (n, m)
    output_dir = os.path.join(script_dir, '../../result/n0n1_grid', tag)
    os.makedirs(output_dir, exist_ok=True)
    for db_name, record_name, record_type in all_db_records:
      ft = pool.submit(task_groundtruth, db_name=db_name,
                       record_name=record_name, record_type=record_type, m=m)
      ft.add_done_callback(callback_groundtruth)
      list_futures.append(ft)
  while True:
    sleep(10)