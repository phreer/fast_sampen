import glob
from scipy.io import savemat
import matplotlib.pyplot as plt

INSTANCE_LINE_NUM = 35 # 每个实验样例 35 行
OFFSET_DATA_LENGTH = 6 - 1
OFFSET_TIME_KD = 19 - 1
OFFSET_TIME_SAMPLING = 33 - 1
OFFSET_ENTROPY = 17 - 1
OFFSET_AB = 18 - 1
OFFSET_ERRS_ENTROPY = 23 - 1
OFFSET_ERRS_A = 24 - 1
OFFSET_ERRS_B = 25 - 1

input_filename_pattern = 'result/convergence_n_final_m%d_r0.15_211014/convergence_r0.15_m%d_*.txt'


for m in range(4, 7):
  list_entropy = []
  list_avg_entropy = []
  list_a = []
  list_b = []
  list_time_kd = []
  list_time_sampling = []
  list_err_entropy = []
  list_err_a = []
  list_err_b = []
  list_std_entropy = []
  list_data_length = []
  filenames = glob.glob(input_filename_pattern)
  for filename in filenames:
    with open(filename, 'r') as f:
      lines = f.readlines()
    for i in range(len(lines) // INSTANCE_LINE_NUM):
      curr_lines = lines[i * INSTANCE_LINE_NUM: (i + 1) * INSTANCE_LINE_NUM]

      data_length = int(curr_lines[OFFSET_DATA_LENGTH][len('\tdata length: '): -1])
      list_data_length.append(data_length)

      entropy = float(curr_lines[OFFSET_ENTROPY][len('\tentropy: '): -1])
      list_entropy.append(entropy)
      
      time_kd = float(curr_lines[OFFSET_TIME_KD][len('\ttime: '): ])
      list_time_kd.append(time_kd)

      time_sampling = float(curr_lines[OFFSET_TIME_SAMPLING][len('\tmean_computation_time: '): -1])
      list_time_sampling.append(time_sampling)

      errs_entropy = map(curr_lines[OFFSET_ERRS_ENTROPY][len('\terrs_entropy: '): -1].split(','))
      err_entropy = sum(errs_entropy) / len(errs_entropy)
      list_err_entropy.append(err_entropy)
      list_avg_entropy.append(entropy + err_entropy)
      std_entropy = sum(map(lambda x: x ** 2, map(lambda x: x - err_entropy))) / (len(errs_entropy) - 1)
      list_std_entropy.append(std_entropy)
    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(list_entropy)
    ax.plot(list_avg_entropy)
    plt.show()