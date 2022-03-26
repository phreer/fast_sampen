import os
import numpy as np
import pandas as pd
from scipy.io import savemat

inputdir = os.path.join(os.environ['HOME'], 'workspace/entropy/additional_data/uurgeg')
outputdir = inputdir
db_name = 'uurgeg_380'
ranges = ['1951-1960', '1961-1970', '1971-1980', '1981-1990', '1991-2000', '2001-2010', '2011-2020']

def read_temperature_from_file(filename):
  d = pd.read_csv(filename)
  return np.array(d['T'])

if __name__ == '__main__':
  list_t = []
  for r in ranges:
    t = read_temperature_from_file(os.path.join(inputdir, 'uurgeg_380_%s.csv' % r))
    list_t.append(t)
  t = np.concatenate(list_t)
  save_dict = {'uurgeg_t': t}
  savemat(os.path.join(outputdir, '%s.mat' % db_name), save_dict)