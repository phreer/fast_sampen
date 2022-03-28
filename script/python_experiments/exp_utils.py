import os
import sys
from enum import Enum
import numpy as np
import wfdb
from scipy.io import loadmat
import sampen
import sampen2d

class RecordType(Enum):
  MAT = 1
  WFDB = 2
  TABLE = 3

try:
  input_dir = os.environ['SAMPEN_DATA_DIR']
except KeyError as e:
  print('please set environment variable SAMPEN_DATA_DIR')

all_db_records = [
  ('kmni_climate', 'uurgeg_t', RecordType.MAT),
  ('gearbox', 'Miss_30_2_1', RecordType.MAT),
  ('rolling_bearing', 'X110_DE_time', RecordType.MAT),
  ('chfdb', 'chf01', RecordType.WFDB),
  ('ltafdb', '00', RecordType.WFDB),
  ('ltstdb', 's20011', RecordType.WFDB),
  ('mghdb', 'mgh001', RecordType.WFDB),
  ('mit-bih-long-term-ecg-database-1.0.0', '14046', RecordType.WFDB),
  ('pink', 'pink-1m.txt', RecordType.TABLE),
  ('gaussian', 'gaussian-1m.txt', RecordType.TABLE),
  ('uniform', 'uniform-1m.txt', RecordType.TABLE),
  ('chbmit', 'chb07_01.txt', RecordType.TABLE),
  ('RRHealth', 'Health_Filt-time-1583.01.rr_multirecord.txt', RecordType.TABLE),
  ('RRCHF', 'CHF_Filt-time-9643.rr_multirecord.txt', RecordType.TABLE),
  ('RRAF', 'AF_fa002.rr_multirecord.txt', RecordType.TABLE),
]

def read_record(db_name, record_name, record_type, n):
  if record_type == RecordType.MAT:
    d = read_mat_record(db_name, record_name, 0, n)
  elif record_type == RecordType.WFDB:
    d = read_wfdb_record(db_name, record_name, 100000, n)
  else:
    d = read_table_record(db_name, record_name, 1, 0, n)
  return d

def read_wfdb_record(db_name: str, record_name: str, sig_offset: int, n: int):
  record_path = os.path.join(input_dir, db_name, record_name)
  sigs = wfdb.rdsamp(record_path, sig_offset, sig_offset + n)[0][:, 0]
  return sigs

def read_mat_record(db_name: str, record_name: str, sig_offset: int, n: int):
  input_path = os.path.join(input_dir, '%s.mat' % db_name)
  d = loadmat(input_path)
  if not record_name in d.keys():
    raise ValueError('record %s not found in database %s' % (record_name, db_name))
  return d[record_name].flatten()[sig_offset: sig_offset + n]

def read_table_record(db_name: str, record_name: str, channel: int, sig_offset: int, n: int = None):
  data = []
  input_path = os.path.join(input_dir, db_name, record_name)
  with open(input_path) as f:
    lines = f.readlines()
  if n and len(lines) < n + sig_offset:
    raise ValueError('the length of the file (%d) is less than n (%d)' % (len(lines), n))
  for i in range(n):
    p = float(lines[sig_offset + i].split()[channel])
    data.append(p)
  return np.array(data)

def get_statistics(a, real_a):
  m = np.mean(a)
  std = np.std(a)
  m_abs_err = np.mean(np.abs(a))
  d = (real_a + 1e-8)
  return m, std, m_abs_err, m / d, std / d, m_abs_err / d

def estimate_sampen_and_save_statistics(sig, r: float, m: int, sig_offset: int,
                                        n0: int, n1: int, outputdir: str,
                                        output_filename: str,
                                        n_computations: int =50,
                                        groundtruth=None):
  sig = np.array(sig).flatten()
  n = len(sig)
  r_scaled = r * np.sqrt(np.var(sig))
  with open(os.path.join(outputdir, output_filename), 'w') as f:
    if groundtruth is None:
      s = sampen.SampEnRKD(sig, r_scaled, m, 0)
      r_e = s.entropy()
      r_a = s.a_norm()
      r_b = s.b_norm()
    else:
      r_e, r_a, r_b = groundtruth
    print(f'output_filename [{output_filename}]', file=f)
    print(f'n0 [{n0}] n1 [{n1}] m [{m}] r [{r}] n [{n}] line_offset [{sig_offset}]', file=f)
    print(f'[RKD] sampen [{r_e:.6f}] a [{r_a:.6f}] b [{r_b:.6f}]', file=f)
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

def estimate_sampen2d_and_save_statistics(sig, r: float, m: int,
                                          width: int, height: int,
                                          n0: int, n1: int, outputdir: str,
                                          output_filename: str,
                                          n_computations: int =50,
                                          groundtruth=None):
  sig = np.array(sig).flatten()
  n = len(sig)
  r_scaled = r * np.std(sig)
  os.makedirs(outputdir, exist_ok=True)
  with open(os.path.join(outputdir, output_filename), 'w') as f:
    if groundtruth is None:
      s = sampen2d.SampEn2DD(sig, r_scaled, m, width, height, 1, 1, 0)
      r_e = s.entropy()
      r_a = s.a_norm()
      r_b = s.b_norm()
      t = s.time()
    else:
      r_e, r_a, r_b = groundtruth
      t = 0
    print(f'output_filename [{output_filename}]', file=f)
    print(f'n0 [{n0}] n1 [{n1}] m [{m}] r [{r}] n [{n}]', file=f)
    print(f'[Direct] sampen2d [{r_e:.6f}] a [{r_a:.6f}] b [{r_b:.6f}] t [{t:.4e}]', file=f)
    list_err_e = np.zeros(n_computations)
    list_err_a = np.zeros(n_computations)
    list_err_b = np.zeros(n_computations)
    list_time = np.zeros(n_computations)
    for i in range(n_computations):
      s = sampen2d.SampEn2DSamplingD(sig, r_scaled, m, width, height, 1, 1,
                                     n0, n1, r_e, r_a, r_b, True, 0)
      err_e = s.entropy() - r_e
      err_a = s.a_norm() - r_a
      err_b = s.b_norm() - r_b
      list_err_e[i] = err_e
      list_err_a[i] = err_a
      list_err_b[i] = err_b
      list_time[i] = s.time()
      print(f'[SWR] sampen2d [{s.entropy():.6f}] a [{s.a_norm():.6f}] b [{s.b_norm():.6f}] t [{s.time():.4e}]', file=f)
      print(f'[SWR] err_sampen2d [{err_e :.4e}] err_a [{err_a:.4e}] err_b [{err_b:.4e}]', file=f)
    m_e, std_e, m_abs_err_e, r_m_e, r_std_e, r_m_abs_err_e = get_statistics(list_err_e, r_e)
    m_a, std_a, m_abs_err_a, r_m_a, r_std_a, r_m_abs_err_a = get_statistics(list_err_a, r_a)
    m_b, std_b, m_abs_err_b, r_m_b, r_std_b, r_m_abs_err_b = get_statistics(list_err_b, r_b)
    print(f'[SWR] m_e [{m_e:.4e}] std_e [{std_e:.4e}] m_abs_err_e [{m_abs_err_e:.4e}] m_t [{np.mean(list_time):.4e}]', file=f)    
    print(f'[SWR] r_m_e [{r_m_e:.4e}] r_std_e [{r_std_e:.4e}] r_m_abs_err_e [{r_m_abs_err_e:.4e}]', file=f)    
    print(f'[SWR] m_a [{m_a:.4e}] std_a [{std_a:.4e}] m_abs_err_a [{m_abs_err_a:.4e}]', file=f)    
    print(f'[SWR] r_m_a [{r_m_a:.4e}] r_std_a [{r_std_a:.4e}] r_m_abs_err_a [{r_m_abs_err_a:.4e}]', file=f)    
    print(f'[SWR] m_b [{m_b:.4e}] std_b [{std_b:.4e}] m_abs_err_b [{m_abs_err_b:.4e}]', file=f)    
    print(f'[SWR] r_m_b [{r_m_b:.4e}] r_std_b [{r_std_b:.4e}] r_m_abs_err_b [{r_m_abs_err_b:.4e}]', file=f)    

# Test
if __name__ == '__main__':
  for db, record, t in all_db_records:
    d = read_record(db, record, t, 10)
    print(db, record, d)