import os
import re

from matplotlib import rcParams
rcParams['font.family'] = 'serif'
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

NUM_LINE_START = 1 
NUM_LINE_CASE = 59
LINE_DATA_LENGTH = 5
LINE_TIME_KD = 20 - 1 - NUM_LINE_START
LINE_TIME_UNIFORM = LINE_TIME_KD + 6
LINE_TIME_SWR = LINE_TIME_UNIFORM + 8
LINE_TIME_QMC = LINE_TIME_SWR + 8
LINE_TIME_QMC_PRESORT = LINE_TIME_QMC + 8
LINE_TIME_FAST_DIRECT = LINE_TIME_QMC_PRESORT + 8
OFFSET_TIME = len('\ttime: ')
ERR_TEMPLATE = '\terror: (-?\d+\.\d+e[+-]\d+|-?nan|-?inf)\terror \(relative\): (-?\d+\.\d+e[+-]\d+|-?nan|-?inf)'
err_pattern = re.compile(ERR_TEMPLATE)

r = 0.1
m = 5
inputdir = 'result/grid_m3_r0.1_210112/time_n04000_n140_r{:.1f}-m{}'.format(r, m)
filename = 'l_100000-{}.txt_2021-01-18.txt'
records = ['00', '14046', 'chf01', 'mgh001', 's20011'] # , 'gaussian_noise-2000000', 'pink_noise-2000000']
outputdir = os.path.join(inputdir, 'fig', 'time_linear')
os.makedirs(outputdir, exist_ok=True)
dpi = 60
figsize=[8, 8]

for record in records:
    input_filename = os.path.join(inputdir, filename.format(record))
    with open(input_filename) as f:
        lines = f.readlines()
    
    lines = lines[NUM_LINE_START: ]
    results = {}
    instances = ['kd', 'uniform', 'swr', 'qmc', 'qmc_presort', 'fast_direct']
    for instance in instances:
        results[instance] = {
            "err_sampen": list(),
            "time": list(), 
            'n': list(), 
        }
    data_length = list()
    i = 0
    while (True):
        curr = lines[NUM_LINE_CASE * i: NUM_LINE_CASE * (i + 1)]
        i = i + 1
        if len(curr) != NUM_LINE_CASE:
            break
        n = int(curr[LINE_DATA_LENGTH][len('\tdata length: '): -1])
        time_kd = float(curr[LINE_TIME_KD][OFFSET_TIME: -1])
        time_uniform = float(curr[LINE_TIME_UNIFORM][OFFSET_TIME: -1])
        err_uniform = float(err_pattern.match(curr[LINE_TIME_UNIFORM + 1]).group(2))
        time_swr = float(curr[LINE_TIME_SWR][OFFSET_TIME: -1])
        err_swr = float(err_pattern.match(curr[LINE_TIME_SWR + 1]).group(2))
        time_qmc = float(curr[LINE_TIME_QMC][OFFSET_TIME: -1])
        err_qmc = float(err_pattern.match(curr[LINE_TIME_QMC + 1]).group(2))
        time_qmc_presort = float(curr[LINE_TIME_QMC_PRESORT][OFFSET_TIME: -1])
        err_qmc_presort = float(err_pattern.match(curr[LINE_TIME_QMC_PRESORT + 1]).group(2))
        time_fast_direct = float(curr[LINE_TIME_FAST_DIRECT][OFFSET_TIME: -1])
        
        data_length.append(n)
        results['kd']['n'].append(n)
        results['kd']['time'].append(time_kd)
        results['kd']['err_sampen'].append(0)
        results['uniform']['n'].append(n)
        results['uniform']['time'].append(time_uniform)
        results['uniform']['err_sampen'].append(err_uniform)
        results['swr']['n'].append(n)
        results['swr']['time'].append(time_swr)
        results['swr']['err_sampen'].append(err_swr)
        results['qmc']['n'].append(n)
        results['qmc']['time'].append(time_qmc)
        results['qmc']['err_sampen'].append(err_qmc)
        results['qmc_presort']['n'].append(n)
        results['qmc_presort']['time'].append(time_qmc_presort)
        results['qmc_presort']['err_sampen'].append(err_qmc_presort)
        results['fast_direct']['n'].append(n)
        results['fast_direct']['time'].append(time_fast_direct)
        results['fast_direct']['err_sampen'].append(0)
    
    print('Parsing Done.')
    savefig_options = {'bbox_inches': 'tight'}
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)
    for instance, result in results.items():
        ax.plot(data_length, result['time'], '-x', label=instance)
        ax.set_ylabel('Time (Seconds)')
        ax.set_yscale('log')
        ax.set_xlabel('Data Length')
        ax.set_xscale('log')
        ax.legend()
    fig.savefig(os.path.join(outputdir, 'time_%s.pdf' % (record)), **savefig_options)
    plt.show()
    # plt.close(fig)