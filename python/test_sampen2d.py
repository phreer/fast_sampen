import os
import sys
import sampen2d
import numpy as np
from PIL import Image, ImageOps

r = 0.3
m = 1
sample_rate = 0.05
n1 = 20

if __name__ == '__main__':
  input_dir = os.path.join(os.environ['HOME'], 'workspace/entropy/2d')
  image_path = os.path.join(input_dir, 'Original_Brodatz/D1.gif')
  image = Image.open(image_path)
  image = ImageOps.grayscale(image)
  image = np.array(image)
  r_scaled = r * np.std(image)
  h, w = image.shape
  n = h * w

  # s = sampen2d.SampEn2DD(image.flatten(), r_scaled, m, w, h, 1, 1, 0)
  # print('time: ', s.time())
  # print('entropy: ', s.entropy())
  # r_e, r_a, r_b = s.entropy(), s.a_norm(), s.b_norm()
  r_e, r_a, r_b = 0, 0, 0

  n0 = n * sample_rate
  s = sampen2d.SampEn2DSamplingD(image.flatten(), r_scaled, m, w, h, 1, 1,
                                 n0, n1, r_e, r_a, r_b, True, 0)
  print('time: ', s.time())
  print('entropy: ', s.entropy())

  s = sampen2d.SampEn2DSamplingD(image.flatten(), r_scaled, m, w, h, 1, 1,
                                 n0, n1, r_e, r_a, r_b, True, 0)
  print('time: ', s.time())
  print('entropy: ', s.entropy())