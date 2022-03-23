import numpy as np
import sampen
import sampen2d

if __name__ == '__main__':
  n = 100000
  m = 2
  r = 0.15
  d = np.random.randn(n)
  r = r * np.sqrt(np.var(d))
  s = sampen.SampEnRKD(d, r, m, 0)
  print('SampEn: %.4f' % s.entropy())
  print('Time: %.4f' % s.time())
  
  n0 = 2000
  n1 = 20
  s = sampen.SampEnSamplingD(d, r, m, n0, n1,
                            s.entropy(), s.a_norm(), s.b_norm(),
                            sampen.RandomType.SWR_UNIFORM, False, False, 0)
  print('SampEn: %.4f' % s.entropy())
  print('Time: %.4f' % s.time())
  
  h, w = 200, 200
  image = np.random.randn(h * w)
  print('Testing SampEn2DD...')
  s = sampen2d.SampEn2DD(image, r, m, w, h, 1, 1, 0)
  print('SampEn: %.4f' % s.entropy())
  print('Time: %.4f' % s.time())

  print('Testing SampEn2DSamplingD...')
  s = sampen2d.SampEn2DSamplingD(image, r, m, w, h, 1, 1, n0, n1, s.entropy(), s.a_norm(), s.b_norm(), 0)
  print('SampEn: %.4f' % s.entropy())
  print('Time: %.4f' % s.time())