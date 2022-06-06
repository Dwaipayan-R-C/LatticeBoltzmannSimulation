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
    
    def visualize_couette(i):
        j= 0
        """Visual function for animate [Refer to animate]"""
        axes[1].plot(couette_velocity_list[i], np.arange(Ny), color = 'r')
    
    def animate(velocity):
        """Creates density animation"""
        anim = animation.FuncAnimation(figs[1],visualize_couette,repeat=True,frames=len(couette_velocity_list))     
        
        anim.save('couette_animation.gif',writer='imagemagic', fps=2)

    """ Calculates the shear wave for Nx by Ny D2Q9 lattice"""
    plt.rcParams.update({'font.size': 16})
    plt.rcParams["font.family"] = "Times New Roman"

    fig1, ax1 = plt.subplots() # Assigned to continuous flow of velocity  
    fig2, ax2 = plt.subplots() # Assigned to continuous flow of velocity  
    figs, axes = [fig1, fig2], [ax1,ax2] 

    
    # Variable declaration
    common_path = os.path.join(output_dir, 'Couette_flow')    
    os.makedirs(common_path, exist_ok=True)
    couette_velocity_list = []

    """Starts from here"""
    size_x = 50                     # 50
    size_y = 50     #+ topbottom_boundary   # 52
      # initilization of the grids used
    f = np.ones((Ny,Nx,9))
    density = np.zeros((Ny,Nx))
    velocity = np.zeros((Ny,Nx,2))
    

    # f = np.ones((9,size_x,size_y))
    # rho_v = np.zeros((size_x,size_y))
    # ux_v = np.zeros((size_x,size_y))
    # uy_v = np.zeros((size_x,size_y))
        # Initializes density and velocity
    # density = np.zeros((Nx,Ny))
    # velocity = np.zeros((2, Nx,Ny))

    # initialize density function w.r.t. rho and u 
    # f = lbm.calculate_equilibrium(density,velocity)    
    # f = np.ones((9,Nx,Ny))

    for step in range(steps):
        print(f'{step+1}//{steps}', end="\r")     
        # print(velocity[:,Nx//2,1])
        f = lbm.stream(f)
        f = lbm.bounce_back(f, lid_vel)        
        rho_v, ux,uy, f = lbm.calculate_collision(f, omega)
        # ux,uy = velocity[:,:,1], velocity[:,:,0]
        
        
        #        
        
            
        # vel_max = max(velocity[:,Nx//2,1])
        # vel_min = min(velocity[:,Nx//2,1])
        
        # # plt.plot(velocity[:,Nx//2,1], np.arange(Ny))
        # # plt.show()
        # # print('')

        if (step%save_every == 0 and step !=0):
            # print('')
            # print(np.min(ux[25,:]))
            # print('')
        #     print(velocity[:,Nx//2,1])
            # axes[0].cla()
            # axes[0].set_ylabel("Width")
            # axes[0].set_xlabel("velocity in X direction")
            # axes[0].plot(ux[25,:], np.arange(Ny))
            # save_path = os.path.join(common_path, f'velocity_at{step}.png')
            # axes[0].set_title('Velocity in x direction for lid velocity {} after {} iteration'.format(lid_vel, step))
            # figs[0].savefig(save_path, bbox_inches='tight', pad_inches=0)
            couette_velocity_list.append(uy[:,25])
        #     print(velocity[:,Nx//2,1])

    animate(couette_velocity_list)
    
          
if __name__ == "__main__":
    couette_flow_simulation(50,50,0.5,'results', 500,5000,5)

