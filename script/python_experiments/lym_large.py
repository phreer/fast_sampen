import sys
import os

import numpy as np
from PIL import Image, ImageOps
import sampen2d

Image.MAX_IMAGE_PIXELS = None

r = 0.3
m = 1
n1 = 20
n0 = 10000

if __name__ == '__main__':
  script_dir = os.path.dirname(sys.argv[0])
  image_path = os.path.join(script_dir, '../../../entropy/data/N55442_5C，6-2021-07-03_03_03_58.tif')
  image = Image.open(image_path)
  image = ImageOps.grayscale(image)
  image = np.array(image)
  r_scaled = r * np.std(image)
  h, w = image.shape
  n = h * w
  s = sampen2d.SampEn2DSamplingD(image.flatten(), r_scaled, m, w, h, 1, 1,
                                 n0, n1, r_e, r_a, r_b, 0)
  print('time: ', s.time())
  print('entropy: ', s.entropy())