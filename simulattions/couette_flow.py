"""Library Imports"""
from turtle import color
import matplotlib.pyplot as plt
import numpy as np
from lbm_common import lbm 
from lbm_common import constant as CV
import os
import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt
from matplotlib import animation
from lbm_common import boundary as boundary

def couette_flow_simulation(Nx: int, Ny: int, omega: float, output_dir: str, save_every, steps, lid_vel):
    """ Calculates the couette flow for Nx by Ny D2Q9 lattice"""

    def visualize_couette(i):
        j= 0
        """Visual function for animate [Refer to animate]"""
        axes[1].plot(couette_velocity_list[i], np.arange(Ny), color = 'darkgreen')
        if(j == i):
            xmin, xmax, ymin, ymax = axes[1].axis()
            x = np.linspace(xmin, xmax)
            y_max = np.full(len(x),ymax-2.5)
            y_min = np.full(len(x),ymin+2.5)
            axes[1].set_ylabel("Width Ny")
            axes[1].set_xlabel("Velocity in X direction")
            axes[1].set_title('Couette Flow with lid velocity {}'.format(lid_vel))
            # setup for gif      
            axes[1].plot(x, y_min, color='k', linewidth=6.0)
            axes[1].plot(x, y_max, color='r', linewidth=6.0)
            axes[1].legend(['Analytical Flow','Simulated Flow','Rigid wall','Moving wall'])
        else:
            j = j+1
    
    
    def analytical_couette(y):
        """Define Analytical function for couette flow"""
        x_vel = lid_vel * y/Ny        
        return x_vel


    def animate(velocity):
        """Creates density animation"""
        anim = animation.FuncAnimation(figs[1],visualize_couette,repeat=True,frames=len(couette_velocity_list))     
        anim.save('{}/Couette_animation.gif'.format(output_dir),writer='imagemagic', fps=2)

    """ Calculates the couette flow for Nx by Ny D2Q9 lattice"""
    plt.rcParams.update({'font.size': CV.fontsize})
    plt.rcParams["font.family"] = CV.fontfamily

    fig1, ax1 = plt.subplots() # Used to save every frame  
    fig2, ax2 = plt.subplots() # Used to create a gif
    figs, axes = [fig1, fig2], [ax1,ax2]
    
    # Variable declaration
    common_path = os.path.join(output_dir, 'Couette_flow')    
    os.makedirs(common_path, exist_ok=True)
    couette_velocity_list = []

    """Starts from here"""
    # initilization of the grids used    
    f = np.ones((Ny,Nx,9))
    velocity = np.zeros((Ny,Nx,2))
    analytical_vel = analytical_couette(np.arange(Ny))
    axes[0].plot(analytical_vel,np.arange(Ny), color = 'purple', linestyle='dashed')
    axes[1].plot(analytical_vel,np.arange(Ny), color = 'purple', linestyle='dashed')
        
    # Iteration starts from here
    for step in range(steps):
        print(f'{step+1}//{steps}', end="\r")  

        # Streaming, Bounceback and Collision
        f = lbm.streaming(f)
        f = boundary.couette_bounce_back(f,lid_vel,velocity)        
               
        f, density, velocity = lbm.calculate_collision(f, omega)        
        # f, density, velocity = lbm.calculate_collision(f, omega)        

        # Saving steps 
        if save_every is not None and (not (step % save_every) and step!=0):
            axes[0].cla()
            axes[0].set_ylabel("Width")
            axes[0].set_xlabel("velocity in X direction")
            axes[0].plot(velocity[:,Nx//2,1], np.arange(Ny), color = 'r')
            save_path = os.path.join(common_path, f'velocity_at{step}.png')
            axes[0].set_title('Velocity in x direction for lid velocity {} after {} iteration'.format(lid_vel, step))
            figs[0].savefig(save_path, bbox_inches='tight', pad_inches=0)
            axes[0].legend(['Analytical','Simulated'])
            couette_velocity_list.append(velocity[:,Nx//2,1])        

    # Animate the couette flow 
    animate(couette_velocity_list)


