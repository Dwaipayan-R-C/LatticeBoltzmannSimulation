"""Library Imports"""
import matplotlib.pyplot as plt
import numpy as np
import lbm 
import constant as CV
import os
from typing import Optional
import sys
import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import make_interp_spline
from matplotlib import animation
from scipy.ndimage.filters import gaussian_filter1d

def couette_flow_simulation(Nx: int, Ny: int, omega: float, output_dir: str, save_every, steps, lid_vel):
    """ Calculates the shear wave for Nx by Ny D2Q9 lattice"""
    plt.rcParams.update({'font.size': 16})
    plt.rcParams["font.family"] = "Times New Roman"

    fig1, ax1 = plt.subplots() #Assigned to continuous flow of velocity  
    figs, axes = [fig1], [ax1] 

    # Variable declaration
    common_path = os.path.join(output_dir, 'Couette_flow')    
    os.makedirs(common_path, exist_ok=True)

    # Initializes density and velocity
    density = np.ones((Ny, Nx))
    velocity_field = np.zeros((Ny, Nx, 2))

    # initialize density function w.r.t. rho and u 
    f = lbm.calculate_equilibrium(density, velocity_field)
    

    for step in range(steps):
        print(f'{step+1}//{steps}', end="\r")      
        
        f = lbm.streaming(f,Ny, Nx, lid_vel)
        f, density, velocity = lbm.collision(f, omega)

        if save_every is not None and (not (step % save_every) ):
            axes[0].cla()
            axes[0].set_ylabel("Width")
            axes[0].set_xlabel("velocity in X direction")
            axes[0].plot(velocity[:,Nx//2,1], np.arange(Ny))
            save_path = os.path.join(common_path, f'velocity_at{step}.png')
            axes[0].set_title('Velocity in x direction for lid velocity {} after {} iteration'.format(lid_vel, step))
            figs[0].savefig(save_path, bbox_inches='tight', pad_inches=0)
          


