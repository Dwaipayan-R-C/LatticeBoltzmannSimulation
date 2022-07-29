# imports
from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
from lbm_common import lbm
from lbm_common import boundary
import os
from matplotlib import animation


def poiseuille_simulation(rho_null, p_diff, output_dir, Nx, Ny, relaxation,steps, save_every):
    """ Calculates the poiseuille simulation for Nx by Ny D2Q9 lattice"""

    def visualize_poiseuille(i):
        j= len(poiseuille_list)
        """Visual function for animate [Refer to animate]"""
        axes[1].plot(poiseuille_list[i], np.arange(Ny))
        
    def analytical_poiseuille():
        x = np.arange(0, Nx+2)
        y = np.arange(0, Ny+2)
        X, Y = np.meshgrid(x, y)
        delta = 2.0 * p_diff /Nx / shear_viscosity / 2.
        y = np.linspace(0, Ny, Ny)
        u_analytical = delta * y * (Ny - y) / 3.
        return u_analytical

    def animate(velocity, vel_max):
        """Creates density animation"""
        axes[1].set_ylabel("Width Ny")
        axes[1].set_xlabel("Velocity in X direction")
        axes[1].set_title('Poiseuille Flow with pressure gradient {}'.format(p_diff))        
        analytical_value = analytical_poiseuille() 
        axes[1].legend(['Analytical Flow','Simulated Flow','Rigid wall','Rigid wall'])
        x = np.linspace(0, np.max(analytical_value), Ny)    
        y_max = np.full(len(x),Ny-1)
        y_min = np.full(len(x),0)   
        axes[1].plot(x, y_max, color='r', linewidth=1.5)
        axes[1].plot(x, y_min, color='r', linewidth=1.5)               
        axes[1].plot(analytical_value, np.arange(len(analytical_value)),color='black', linestyle='dashed')
        anim = animation.FuncAnimation(figs[1],visualize_poiseuille,repeat=True,frames=len(poiseuille_list))        
        anim.save('{}/Poiseuille_animation.gif'.format(output_dir),writer='imagemagic', fps=2)

    # Variable declaration
    common_path = os.path.join(output_dir, 'Poiseuille_flow')    
    os.makedirs(common_path, exist_ok=True)
    plt.rcParams.update({'font.size': 10})
    plt.rcParams["font.family"] = "Cambria"

    fig1, ax1 = plt.subplots() # Used to save every frame  
    fig2, ax2 = plt.subplots() # Used to create a gif
    fig3, ax3 = plt.subplots() # Used to create a gif
    figs, axes = [fig1, fig2,fig3], [ax1,ax2,ax3]
    shear_viscosity = (1/relaxation-0.5)/3
    
    poiseuille_list = []        
    rho_in = rho_null+p_diff
    rho_out = rho_null-p_diff
    velocity = np.zeros((Ny+2,Nx+2,2))
    # initialize
    rho = np.ones((Ny+2, Nx + 2))    
    f = lbm.calculate_equilibrium(rho,velocity)

    # loop
    for step in range(steps):
        print(f'{step+1}//{steps}', end="\r")
        rho, velocity, f =lbm.periodic_boundary_with_pressure_variations(f,rho_in,rho_out)
        f = lbm.streaming(f)
        # Apply boundary        
        f = boundary.f_rigid_wall(f,True,True,False,False)
        # f = boundary.poiseuille_bounce_back(f,0)        
        velocity = lbm.calculate_velocity(f,rho)
        f, density, velocity = lbm.calculate_collision(f, relaxation)
        if save_every is not None and (not (step % save_every) and step!=0):
             
             
            axes[0].cla()
            axes[0].set_ylabel("Width of the channel")
            axes[0].set_xlabel("Velocity [m/s] in X direction")
            x_val = velocity[1:-1, Nx//2,1]            
            y_val = np.arange(Ny)            
            axes[0].plot(x_val, y_val, color = 'g')           
            
            save_path = os.path.join(common_path, f'velocity_at{step}.png')
            axes[0].set_title('Poiseuille Flow with pressure gradient {}'.format(p_diff))
            analytical_value = analytical_poiseuille()       
        
            x = np.linspace(0, np.max(analytical_value), Ny)    
            y_max = np.full(len(x),Ny-1)
            y_min = np.full(len(x),0)
            axes[0].plot(x, y_max, color='r', label="Fixed wall")
            axes[0].plot(x, y_min, color='r')
            axes[0].plot(analytical_value, np.arange(len(analytical_value)),color='black', linestyle='dashed', label='Analytical')
            axes[0].legend(['Analytical','Simulated'])
            figs[0].savefig(save_path,  pad_inches=1)            
            poiseuille_list.append(velocity[1:-1, Nx//2,1])
        

    # visualize
    analytical_poiseuille()
    animate(velocity,np.max(velocity[1:-1, Nx//2,1]))
    x, y = np.meshgrid(np.arange(Nx+2), np.arange(Ny+2))
    save_density_path = os.path.join(common_path, f'density.png')
    ax3.scatter(x,y,c=density)
    figs[2].savefig(save_density_path,  pad_inches=1)