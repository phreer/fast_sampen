from genericpath import exists
import os
import numpy as np
import sampen
import sampen2d

def read_record(input: str, input_format: str, n: int):
  data = []
  with open(input) as f:
    lines = f.readlines()
  if n and len(lines) <= n:
    raise ValueError('the length of the file (%d) is less than n (%d)' % (len(lines), n))
  for i in range(n):
    if input_format == 'simple':
      p = float(lines[i])
    else:
     p = float(lines[i].split()[1])
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
    else:
      r_e, r_a, r_b = groundtruth
    print(f'output_filename [{output_filename}]', file=f)
    print(f'n0 [{n0}] n1 [{n1}] m [{m}] r [{r}] n [{n}]', file=f)
    print(f'[Direct] sampen2d [{r_e:.6f}] a [{r_a:.6f}] b [{r_b:.6f}]', file=f)
    list_err_e = np.zeros(n_computations)
    list_err_a = np.zeros(n_computations)
    list_err_b = np.zeros(n_computations)
    for i in range(n_computations):
      s = sampen2d.SampEn2DSamplingD(sig, r_scaled, m, width, height, 1, 1,
                                     n0, n1, r_e, r_a, r_b, 0)
      err_e = s.entropy() - r_e
      err_a = s.a_norm() - r_a
      err_b = s.b_norm() - r_b
      list_err_e[i] = err_e
      list_err_a[i] = err_a
      list_err_b[i] = err_b
      print(f'[SWR] sampen2d [{s.entropy():.6f}] a [{s.a_norm():.6f}] b [{s.b_norm():.6f}]', file=f)
      print(f'[SWR] err_sampen2d [{err_e :.4e}] err_a [{err_a:.4e}] err_b [{err_b:.4e}]', file=f)
    m_e, std_e, m_abs_err_e, r_m_e, r_std_e, r_m_abs_err_e = get_statistics(list_err_e, r_e)
    m_a, std_a, m_abs_err_a, r_m_a, r_std_a, r_m_abs_err_a = get_statistics(list_err_a, r_a)
    m_b, std_b, m_abs_err_b, r_m_b, r_std_b, r_m_abs_err_b = get_statistics(list_err_b, r_b)
    print(f'[SWR] m_e [{m_e:.4e}] std_e [{std_e:.4e}] m_abs_err_e [{m_abs_err_e:.4e}]', file=f)    
    print(f'[SWR] r_m_e [{r_m_e:.4e}] r_std_e [{r_std_e:.4e}] r_m_abs_err_e [{r_m_abs_err_e:.4e}]', file=f)    
    print(f'[SWR] m_a [{m_a:.4e}] std_a [{std_a:.4e}] m_abs_err_a [{m_abs_err_a:.4e}]', file=f)    
    print(f'[SWR] r_m_a [{r_m_a:.4e}] r_std_a [{r_std_a:.4e}] r_m_abs_err_a [{r_m_abs_err_a:.4e}]', file=f)    
    print(f'[SWR] m_b [{m_b:.4e}] std_b [{std_b:.4e}] m_abs_err_b [{m_abs_err_b:.4e}]', file=f)    
    print(f'[SWR] r_m_b [{r_m_b:.4e}] r_std_b [{r_std_b:.4e}] r_m_abs_err_b [{r_m_abs_err_b:.4e}]', file=f)    