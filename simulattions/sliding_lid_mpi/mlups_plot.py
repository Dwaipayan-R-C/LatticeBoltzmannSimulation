
"""Imports"""
import matplotlib.pyplot as plt
import numpy as np
import os

"""Code runs"""
# plt info and directory management
data_files = ['output_200.out','output_400.out','output_600.out','output_800.out','output_1000.out']
cuurent_dir = os.getcwd()
save_path = os.path.join(cuurent_dir,'simulattions\\sliding_lid_mpi\\data')
input_file = os.path.join(save_path,)
plt.rcParams["figure.figsize"] = (10,5.5)
fig,ax = plt.subplots()

# Runs through all the files
for file in data_files:
    with open(os.path.join(save_path,f'{file}'))as f:
        lines = f.readlines()
        
    list_lines = []
    for i in range(1,8):
        split_list = lines[i].split(' ')
        updated_string = []
        for j in split_list:
            if(j!=''):
                j.splitlines()
                updated_string.append(int(j))
    
        list_lines.append(updated_string)

    grid_array = np.array(list_lines)
    ml_ups = grid_array[:,-1]
    seconds = grid_array[:,-2]
    grid_len = grid_array[:,0]
    processes = grid_array[:,1]**2
    ax.plot(processes,ml_ups, '-o')
    
# Plot structuring
ax.legend(['200 x 200','400 x 400','600 x 600','800 x 800','1000 x 1000'])
ax.set_title('MLups (seconds) vs No. of processes for different grid size')
ax.grid()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('No. of processes')
ax.set_ylabel('MLUps in seconds')
save_path = os.path.join(f'results/MLups_plot.png') 
fig.savefig(save_path, pad_inches=100)

