import sys
import os
from os.path import join as pjoin
from concurrent.futures import ProcessPoolExecutor


import numpy as np
from PIL import Image, ImageOps
import sampen2d
from exp_utils import estimate_sampen2d_and_save_statistics

Image.MAX_IMAGE_PIXELS = None

n1 = 20
n0 = 10000
n_workers = 35

list_db = ['Original_Brodatz', 'Colored_Brodatz/', 'Normalized_Brodatz',
           'Multi_band_Texture_database/', 'original_image']
input_dir = pjoin(os.environ['HOME'], 'workspace/entropy/2d')
script_dir = os.path.dirname(sys.argv[0])

def estimate_sampen2d_and_save_task(db_name, record_name, r, m):
  output_dir = pjoin(script_dir, '../../result/2d_err/m%d_r%.2f' % (m, r), db_name)
  output_filename = os.path.splitext(record_name)[0] + '.txt'
  input_path = pjoin(input_dir, db_name, record_name)
  image = Image.open(input_path)
  image = np.array(ImageOps.grayscale(image))
  height, width = image.shape
  print('start %s %s %.2f %d...' % (db_name, record_name, r, m))
  estimate_sampen2d_and_save_statistics(image.flatten(), r, m,
                                        width=width, height=height,
                                        n0=n0, n1=n1, outputdir=output_dir,
                                        output_filename=output_filename)
  
if __name__ == '__main__':
  pool = ProcessPoolExecutor(max_workers=n_workers)
  for r in [0.15, 0.3]:
    for m in [1, 2, 3]:
      for db in list_db:
        for record in os.listdir(pjoin(input_dir, db)):
          pool.submit(estimate_sampen2d_and_save_task, db_name=db, record_name=record, r=r, m=m)
