"""Library Imports"""

import matplotlib.pyplot as plt
# import ipyparallel as ipp
import numpy as np
from lbm_common import lbm 
from lbm_common import constant as CV
import os
import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt
from matplotlib import animation
from lbm_common import boundary as boundary

# def mpi_example():
#     from mpi4py import MPI
#     comm = MPI.COMM_WORLD
#     return f"Hello World from rank {comm.Get_rank()}. total ranks={comm.Get_size()}"

# request an MPI cluster with 4 engines


def sliding_lid_simulation(Nx: int, Ny: int, re: float, output_dir: str, save_every, steps, lid_vel):
    """ Calculates the sliding_lid flow for Nx by Ny D2Q9 lattice"""

    omega = (2*re)/(6*Nx*lid_vel+re)
    def visualize_sliding_lid(i):
        
        """Visual function for animate [Refer to animate]"""
        speed = np.sqrt(velocity[:, :, 0].T**2 + velocity[:, :, 1].T**2)
        x, y = np.meshgrid(np.arange(Nx+2), np.arange(Ny+2))
        axes[1].cla()
        axes[1].set_ylabel("Width Y")
        axes[1].set_xlabel("Length X")
        axes[1].set_title('Sliding lid plot for lid velocity {} and reynold number {}'.format(lid_vel, re))
        axes[1].streamplot(x, y, sliding_lid_velocity_list[i][:, :, 1], sliding_lid_velocity_list[i][:, :, 0],color=speed, cmap = plt.cm.jet)         
  

    def animate(velocity):
        """Creates density animation"""
        anim = animation.FuncAnimation(figs[1],visualize_sliding_lid,repeat=True,frames=len(sliding_lid_velocity_list), cache_frame_data = False)     
        anim.save('{}/Karman_vortex_animation.gif'.format(output_dir),writer='imagemagic', fps=2)
        
    # comm = MPI.COMM_WORLD
    # size = comm.Get_size()
    # rank = comm.Get_rank()

    # CommCart = comm.Create_cart((1, 1), periods=(False, False))
    
    # print("The rank is: {} and the coordinate is: {} ".format(rank, CommCart.Get_coords(rank)))
    
    """ Calculates the sliding_lid flow for Nx by Ny D2Q9 lattice"""
    plt.rcParams.update({'font.size': CV.fontsize})
    plt.rcParams["font.family"] = CV.fontfamily

    fig1, ax1 = plt.subplots() # Used to save every frame  
    fig2, ax2 = plt.subplots() # Used to create a gif
    figs, axes = [fig1, fig2], [ax1,ax2]
    
    # Variable declaration
    common_path = os.path.join(output_dir, 'Sliding_lid')    
    os.makedirs(common_path, exist_ok=True)
    sliding_lid_velocity_list = []

    """Starts from here"""
    # initilization of the grids used    
    density = np.ones((Ny+2,Nx+2))
    velocity = np.zeros((Ny+2,Nx+2,2))
    f = lbm.calculate_equilibrium(density,velocity)
    
        
    # Iteration starts from here
    for step in range(steps):
        print(f'{step+1}//{steps}', end="\r")  

        # Streaming, Bounceback and Collision
        f = lbm.streaming(f)
        f = boundary.sliding_bounce_back(f,lid_vel,velocity)
        f, density, velocity = lbm.calculate_collision(f, omega)        
        
        # Saving steps 
        if save_every is not None and (not (step % save_every)) and step !=0:
            axes[0].cla()
            axes[0].set_ylabel("Width of grid (Y)")
            axes[0].set_xlabel("Length of grid (X)")
            x, y = np.meshgrid(np.arange(Nx+2), np.arange(Ny+2)) 
            speed = np.sqrt(velocity[:, :, 0].T**2 + velocity[:, :, 1].T**2)           
            axes[0].streamplot(x, y, velocity[:, :, 1], velocity[:, :, 0], color=speed, cmap = plt.cm.jet)
            save_path = os.path.join(common_path, f'sliding_{step}.png')
            axes[0].set_title('Sliding lid flow with lid velocity {} and reynolds number {} after {} iteration'.format(lid_vel,re, step))
            figs[0].savefig(save_path, bbox_inches='tight', pad_inches=0)            
            sliding_lid_velocity_list.append(velocity)        

    # # Animate the sliding_lid flow 
    
    animate(sliding_lid_velocity_list)


    
