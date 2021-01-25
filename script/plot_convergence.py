import os

from matplotlib import rcParams
rcParams['font.family'] = 'serif'
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

NUM_LINE_SETTING = 19
NUM_LINE_CASE = 58
LINE_SAMPLE_SIZE = 1
LINE_SAMPLE_NUM = LINE_SAMPLE_SIZE + 1
LINE_ERR_SAMPEN_KD = LINE_SAMPLE_SIZE + 5
LINE_GAP_INSTANCE = 9
LINE_ERR_SAMPEN_SWR = LINE_ERR_SAMPEN_KD + LINE_GAP_INSTANCE
LINE_ERR_SAMPEN_QMC = LINE_ERR_SAMPEN_SWR + LINE_GAP_INSTANCE
LINE_ERR_SAMPEN_QMC_PRESORT = LINE_ERR_SAMPEN_QMC + LINE_GAP_INSTANCE
LINE_ERR_SAMPEN_GRID = LINE_ERR_SAMPEN_QMC_PRESORT + LINE_GAP_INSTANCE
LINE_ERR_SAMPEN_GRID_PRESORT = LINE_ERR_SAMPEN_GRID + LINE_GAP_INSTANCE
LINE_OFFSET_STD = 1
LINE_OFFSET_A = 2
LINE_OFFSET_A_STD = 3
LINE_OFFSET_B = 4
LINE_OFFSET_B_STD = 5

inputdir = 'result/grid_m3_r0.1_210112/'
filename = 'convergence_r0.1_m3_%s.txt_2021-01-12.txt'
records = ['00', 'chf01', 'mgh001', 's20011', 'gaussian_noise-2000000', 'pink_noise-2000000']
outputdir = os.path.join(inputdir, 'fig', 'convergence')
os.makedirs(outputdir, exist_ok=True)

for record in records:
    input_filename = os.path.join(inputdir, filename % record)
    with open(input_filename) as f:
        lines = f.readlines()
    
    increment_ss = 200
    increment_sn = 10
    lines = lines[NUM_LINE_SETTING: ]
    results = {}
    instances = ['uniform', 'swr', 'qmc', 'qmc_presort', 'grid', 'grid_presort']
    num_ss = 20
    num_sn = 25
    num_case = num_ss * num_sn
    for instance in instances:
        results[instance] = {
            "err_sampen": np.zeros([num_ss, num_sn]), 
            "std_sampen": np.zeros([num_ss, num_sn]), 
            "std_a": np.zeros([num_ss, num_sn]), 
            "err_a": np.zeros([num_ss, num_sn]), 
            "err_b": np.zeros([num_ss, num_sn]), 
            "std_b": np.zeros([num_ss, num_sn]), 
        }
    
    for i in range(num_case - 2):
        curr = lines[NUM_LINE_CASE * i: NUM_LINE_CASE * (i + 1)]
        n0 = int(curr[LINE_SAMPLE_SIZE][13: -1])
        n1 = int(curr[LINE_SAMPLE_NUM][12: -1])
        i_n0 = n0 // increment_ss - 1
        i_n1 = n1 // increment_sn - 1
        for j, instance in enumerate(instances):
            LINE_INSTANCE = LINE_ERR_SAMPEN_KD + LINE_GAP_INSTANCE * j
            err_sampen = float(curr[LINE_INSTANCE][len('\tmean_errs_sampen: '): -1])
            std_sampen = float(curr[LINE_INSTANCE + LINE_OFFSET_STD][len('\tstd_errs_sampen: '): -1])
            err_a = float(curr[LINE_INSTANCE + LINE_OFFSET_A][len('\tmean_errs_a: '): -1])
            std_a = float(curr[LINE_INSTANCE + LINE_OFFSET_A_STD][len('\tstd_errs_a: '): -1])
            err_b = float(curr[LINE_INSTANCE + LINE_OFFSET_B][len('\tmean_errs_b: '): -1])
            std_b = float(curr[LINE_INSTANCE + LINE_OFFSET_B_STD][len('\tstd_errs_b: '): -1])
            results[instance]['err_sampen'][i_n0, i_n1] = err_sampen
            results[instance]['std_sampen'][i_n0, i_n1] = std_sampen
            results[instance]['err_a'][i_n0, i_n1] = err_a
            results[instance]['std_a'][i_n0, i_n1] = std_a
            results[instance]['err_b'][i_n0, i_n1] = err_b
            results[instance]['std_b'][i_n0, i_n1] = std_b
    
    print('Parsing Done.')
    
    n0s = np.arange(200, 4001, 200)
    n1s = np.arange(10, 251, 10)
    n0s, n1s = np.meshgrid(n1s, n0s)
    dpi = 60
    figsize=[8, 6]
    savefig_options = {'bbox_inches': 'tight'}
    for instance in instances:
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot(111, projection='3d')
        print(type(ax))
        surf = ax.plot_surface(n1s, n0s, results[instance]['err_sampen'], cmap=cm.coolwarm)
        ax.set_zlabel('Mean Error of Sample Entropy', labelpad=10)
        ax.set_xlabel('Sample Size')
        ax.set_ylabel('The Number of Computations')
        fig.savefig(os.path.join(outputdir, 'mean_err_sampen_%s_%s.pdf' % (record, instance)), **savefig_options)
        plt.close(fig)
        
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(n1s, n0s, np.abs(results[instance]['std_sampen']), cmap=cm.coolwarm)
        ax.set_zlabel('Standard Deviation of Sample Entropy')
        ax.set_xlabel('Sample Size')
        ax.set_ylabel('The Number of Computations')
        fig.savefig(os.path.join(outputdir, 'std_err_sampen_%s_%s.pdf' % (record, instance)), **savefig_options)
        plt.close(fig)
        
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(n1s, n0s, results[instance]['err_a'], cmap=cm.coolwarm)
        ax.set_zlabel('Mean Error of Matching Probability')
        ax.set_xlabel('Sample Size')
        ax.set_ylabel('The Number of Computations')
        fig.savefig(os.path.join(outputdir, 'mean_err_a_%s_%s.pdf' % (record, instance)), **savefig_options)
        plt.close(fig)
        
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(n1s, n0s, np.abs(results[instance]['std_a']), cmap=cm.coolwarm)
        ax.set_zlabel('Standard Deviation of Matching Probability')
        ax.set_xlabel('Sample Size')
        ax.set_ylabel('The Number of Computations')
        fig.savefig(os.path.join(outputdir, 'std_err_a_%s_%s.pdf' % (record, instance)), **savefig_options)
        plt.close(fig)
        
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(n1s, n0s, results[instance]['err_b'], cmap=cm.coolwarm)
        ax.set_zlabel('Mean Error of Matching Probability')
        ax.set_xlabel('Sample Size')
        ax.set_ylabel('The Number of Computations')
        fig.savefig(os.path.join(outputdir, 'mean_err_b_%s_%s.pdf' % (record, instance)), **savefig_options)
        plt.close(fig)
        
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(n1s, n0s, np.abs(results[instance]['std_b']), cmap=cm.coolwarm)
        ax.set_zlabel('Standard Deviation of Matching Probability')
        ax.set_xlabel('Sample Size')
        ax.set_ylabel('The Number of Computations')
        fig.savefig(os.path.join(outputdir, 'std_err_b_%s_%s.pdf' % (record, instance)), **savefig_options)
        plt.close(fig)