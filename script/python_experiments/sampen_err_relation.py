import os
from statistics import mean
import threading
import sys
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor

import sampen
import wfdb

script_dir = os.path.dirname(sys.argv[0])
inputdir = '/home/shared/data.PhysioNet'
outputdir = os.path.join(script_dir, '../../result/sampen_err_relation')
sig_offset = 100000
n = 200000
m = 2
r = 0.15
n0 = 2048
n1 = 20
n_computations = 50
n_workers = 10

def read_record_file(filename):
  with open(filename) as f:
    records = list(f.readlines())
  return records

def read_record(record_name):
  sigs = wfdb.rdsamp(record_name, sig_offset, sig_offset + n)[0]
  return sigs

def get_statistics(a, real_a):
  m = np.mean(a)
  std = np.std(a)
  m_abs_err = np.mean(np.abs(a))
  d = (real_a + 1e-8)
  return m, std, m_abs_err, m / d, std / d, m_abs_err / d

def compute_sampen_and_save(sig, output_filename):
  with open(os.path.join(outputdir, output_filename), 'w') as f:
    r_scaled = r * np.sqrt(np.var(sig))
    s = sampen.SampEnRKD(sig, r_scaled, m, 0)
    r_e = s.entropy()
    r_a = s.a_norm()
    r_b = s.b_norm()
    print(f'thread_id [{threading.current_thread().name}]: {output_filename}')
    print(f'output_filename [{output_filename}]', file=f)
    print(f'n0 [{n0}] n1 [{n1}] m [{m}] r [{r}] n [{n}] line_offset [{sig_offset}]', file=f)
    print(f'[RKD] sampen [{s.entropy():.6f}] a [{s.a_norm():.6f}] b [{s.b_norm():.6f}]', file=f)
    list_err_e = np.zeros(n_computations)
    list_err_a = np.zeros(n_computations)
    list_err_b = np.zeros(n_computations)
    for i in range(n_computations):
      s = sampen.SampEnSamplingD(sig, r_scaled, m, n0, n1, r_e, r_a, r_b,
                                 sampen.SWR_UNIFORM, True, False, 0)
      err_e = s.entropy() - r_e
      err_a = s.a_norm() - r_a
      err_b = s.b_norm() - r_b
      list_err_e[i] = err_e
      list_err_a[i] = err_a
      list_err_b[i] = err_b
      print(f'[SWR] sampen [{s.entropy():.6f}] a [{s.a_norm():.6f}] b [{s.b_norm():.6f}]', file=f)
      print(f'[SWR] err_sampen [{err_e :.4e}] err_a [{err_a:.4e}] err_b [{err_b:.4e}]', file=f)
    m_e, std_e, m_abs_err_e, r_m_e, r_std_e, r_m_abs_err_e = get_statistics(list_err_e, r_e)
    m_a, std_a, m_abs_err_a, r_m_a, r_std_a, r_m_abs_err_a = get_statistics(list_err_a, r_a)
    m_b, std_b, m_abs_err_b, r_m_b, r_std_b, r_m_abs_err_b = get_statistics(list_err_b, r_b)
    print(f'[SWR] m_e [{m_e:.4e}] std_e [{std_e:.4e}] m_abs_err_e [{m_abs_err_e:.4e}]', file=f)    
    print(f'[SWR] r_m_e [{r_m_e:.4e}] r_std_e [{r_std_e:.4e}] r_m_abs_err_e [{r_m_abs_err_e:.4e}]', file=f)    
    print(f'[SWR] m_a [{m_a:.4e}] std_a [{std_a:.4e}] m_abs_err_a [{m_abs_err_a:.4e}]', file=f)    
    print(f'[SWR] r_m_a [{r_m_a:.4e}] r_std_a [{r_std_a:.4e}] r_m_abs_err_a [{r_m_abs_err_a:.4e}]', file=f)    
    print(f'[SWR] m_b [{m_b:.4e}] std_b [{std_b:.4e}] m_abs_err_b [{m_abs_err_b:.4e}]', file=f)    
    print(f'[SWR] r_m_b [{r_m_b:.4e}] r_std_b [{r_std_b:.4e}] r_m_abs_err_b [{r_m_abs_err_b:.4e}]', file=f)    

def process_record(db_name, record_name):
  record_path = os.path.join(inputdir, db_name, record_name)
  try:
    sigs = read_record(record_path)
  except Exception as e:
    print(e, file=sys.stderr)
    print(f'Failed to read signal [{record_name}].')
    return
  n_sigs = sigs.shape[1]
  for sig_i in range(n_sigs):
    output_filename = '%s_%s_%d.txt' % (db_name, record_name, sig_i)
    sig = sigs[:, sig_i]
    compute_sampen_and_save(sig, output_filename)
  
def process_db(db_name, pool=None):
  folder = os.path.join(inputdir, db_name)
  records = read_record_file(os.path.join(folder, 'RECORDS'))
  for record_name in records:
    record_name = record_name.strip()
    pool.submit(process_record, db_name=db_name, record_name=record_name)


if __name__ == '__main__':
  os.makedirs(outputdir, exist_ok=True)
  dbs = ['chfdb', 'ltstdb', 'ltafdb', 'mghdb', 'mit-bih-long-term-ecg-database-1.0.0']
  pool = ProcessPoolExecutor(max_workers=n_workers)
  for db in dbs:
    process_db(db, pool)