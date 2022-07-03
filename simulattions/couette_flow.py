"""Library Imports"""
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
        axes[1].plot(couette_velocity_list[i], couette_velocity_list_y[i])
        # axes[1].plot(couette_velocity_list[i], np.arange(1,Ny))
        if(j == i):
            xmin, xmax, ymin, ymax = axes[1].axis()
            x = np.linspace(xmin, xmax)
            y_max = np.full(len(x),ymax-2)
            y_min = np.full(len(x),ymin+2)
            axes[1].set_ylabel("Width Ny")
            axes[1].set_xlabel("Velocity in X direction")
            axes[1].set_title('Couette Flow with lid velocity {}'.format(lid_vel))
            # setup for gif      
            # axes[1].plot(x, y_min, color='k')
            # axes[1].plot(x, y_max, color='r')
            # axes[1].legend(['Analytical Flow','Simulated Flow','Rigid wall','Moving wall'])
        else:
            j = j+1
    
    
    def analytical_couette(y):
        """Define Analytical function for couette flow"""
        x_vel = lid_vel * (y-1)/(Ny-2)        
        return x_vel


    def animate(velocity, couette_velocity_list_y):
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
    couette_velocity_list, couette_velocity_list_y = [], []


    """Starts from here"""
    # initilization of the grids used    
    f = np.ones((Ny,Nx,9))
    velocity = np.zeros((Ny,Nx,2))
    analytical_vel = analytical_couette(np.arange(1,Ny-1))    
    axes[1].plot(analytical_vel,np.arange(1,Ny-1), color = 'purple', linestyle='dashed')
        
    # Iteration starts from here
    for step in range(steps):
        print(f'{step+1}//{steps}', end="\r")  

        # Streaming, Bounceback and Collision
        f, density, velocity = lbm.calculate_collision(f, omega) 
        f = lbm.streaming(f)
        f = boundary.couette_bounce_back(f,lid_vel,velocity)        
               
               
        # f, density, velocity = lbm.calculate_collision(f, omega)        

        # Saving steps 
        if save_every is not None and (not (step % save_every) and step!=0):

            x_0 = velocity[-2,Nx//2,1]
            x_1 = velocity[-3,Nx//2,1]
            y_0 = Ny - 2
            y_1 = Ny - 3
            y_ex = Ny - 1.5
            x_ex = (y_ex - (y_0-((y_1-y_0)/(x_1-x_0))*x_0))/((y_1-y_0)/(x_1-x_0))
            low_y = 0.5
            low_x1 = velocity[2,Nx//2,1]
            low_x = low_y * low_x1
            plt_velocity = velocity
            # plt_velocity[0,Nx//2,1] = plt_velocity[1,Nx//2,1]
            # plt_velocity[-1,Nx//2,1] = x_ex        
            
            axes[0].cla()
            axes[0].set_ylabel("Width")
            axes[0].set_xlabel("velocity in X direction")
            
            # x_val = plt_velocity[1:-1,Nx//2,1]
            x_val = plt_velocity[1:-1,Nx//2,1]
            # y_val = np.arange(1,Ny-1)
            y_val = np.arange(1,Ny-1)
            axes[0].plot(analytical_vel,np.arange(1,Ny-1), color = 'purple', linestyle='dashed')
            axes[0].plot(x_val, y_val)            
            # axes[0].plot(velocity[1:-1,Nx//2,1], np.arange(1,Ny-1), color = 'r')
            save_path = os.path.join(common_path, f'velocity_at{step}.png')
            axes[0].set_title('Velocity in x direction for lid velocity {} after {} iteration'.format(lid_vel, step))
            figs[0].savefig(save_path, bbox_inches='tight', pad_inches=0)
            axes[0].legend(['Analytical','Simulated'])
            couette_velocity_list.append(x_val)        
            couette_velocity_list_y.append(y_val)        
            # couette_velocity_list.append(velocity[1:-1,Nx//2,1])        

    # Animate the couette flow 
    animate(couette_velocity_list, couette_velocity_list_y)