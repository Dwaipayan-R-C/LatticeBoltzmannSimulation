"""Library Imports"""
import matplotlib.pyplot as plt
import numpy as np
from lbm_common import lbm
from lbm_common import constant as CV
import os
import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
from scipy.interpolate import make_interp_spline
from matplotlib import animation


def shear_wave_simulation(Nx: int, Ny: int, omega: float, eps: float, output_dir: str, save_every, steps, exp_type: str):
    """ Calculates the shear wave for Nx by Ny D2Q9 lattice"""

    plt.rcParams.update({'font.size': CV.fontsize})
    plt.rcParams["font.family"] = CV.fontfamily
    fig1, ax1 = plt.subplots() #Assigned to continuous flow of velocity    
    fig2, ax2 = plt.subplots() #Assigned to view velocity stream or density flow
    fig3, ax3 = plt.subplots() #Assigned to sinusoidal decay over time
    fig4, ax4 = plt.subplots() #Assigned to intial density / velocity
    figs, axes = [fig1, fig2,fig3,fig4], [ax1, ax2,ax3,ax4]  
    
    def instant_theoretical_velocity(v):
        """Define Analytical Velocity value"""
        y = np.exp(-v*(2*np.pi/Ny) ** 2)
        return y
        
    def visualize_density(i):
        """Visual function for animate [Refer to animate]"""
        axes[3].scatter(x,y,c=density_animation[i], vmin=density_animation[i].min(), vmax=density_animation[i].max())
    
    def animate(density_animation):
        """Creates density animation"""
        axes[3].set_xlabel('X meshgrid')
        axes[3].set_xlabel('Y meshgrid')
        axes[3].set_title('Density flow in shear wave decay')
        anim = animation.FuncAnimation(figs[3],visualize_density,repeat=False,frames=len(density_animation))
        anim.save('{}/Density_animation.gif'.format(output_dir),writer='imagemagic', fps=2)
    
    def decay_perturbation(t, viscosity):
        """Creates decay perturbation"""
        size = Ny if exp_type == CV.velocity else Nx
        return eps * np.exp(-viscosity * (2*np.pi/size)**2 * t)

    # Create Meshgrid
    x, y = np.meshgrid(np.arange(Nx), np.arange(Ny))
    plt.rcParams.update({'font.size': CV.fontsize})
    plt.rcParams["font.family"] = CV.fontfamily

    """Initializes the velocity and density"""
    if exp_type == CV.velocity:
        density = np.ones((Ny, Nx), dtype=np.float32)
        velocity = np.zeros((Ny, Nx, 2), dtype=np.float32)
        velocity[:, :, 1] =  eps * np.sin(2*np.pi/Ny*y)
    else:
        density = 1 + eps * np.sin(2*np.pi/Nx*x)
        velocity = np.zeros((Ny, Nx, 2), dtype=np.float32)

    
    """Variable Declarations"""
    f = lbm.calculate_equilibrium(density, velocity)
    density_animation = []
    den_or_vel, x_value, max_min_list, den_or_vel_list, theoretical_velocity =[], [], [],[],[]
    common_path = os.path.join(output_dir, 'shear_decay')
    v_not = velocity
    

    
    for step in range(steps):
        """Loop starts"""
        print(f'{step+1}//{steps}', end="\r")
        f = lbm.streaming(f)
        f, density, velocity = lbm.calculate_collision(f, omega)
        if step%save_every==0:
            density_animation.append(density) #if density flow is to be observed
            path = os.path.join(common_path, f'decay_{exp_type}')
            os.makedirs(path, exist_ok=True)
            
            """Density / Velocity continuous flow"""
            # 1. Check the continuous velocity flow
            if (exp_type == CV.density):
                axes[1].cla()
                axes[1].scatter(x,y,c=density, vmin=density.min(), vmax=density.max())
                figs[1].canvas.draw()
                figs[1].canvas.flush_events()
                figs[1].show()
            else:
                axes[1].cla()
                v = np.sqrt(velocity[:, :, 0]**2 + velocity[:, :, 1]**2)
                axes[1].imshow(v, cmap='RdBu_r', vmin=0, interpolation='spline16')            
                figs[1].canvas.draw()
                figs[1].canvas.flush_events()
                figs[1].show()


            """Density / Velocity Decay plot and save"""
            if(exp_type == CV.density):          
               
                # 2. Save the sin flow of density at a point            
                axes[0].cla()
                axes[0].set_ylim([-eps + 1, eps + 1])
                axes[0].plot(np.arange(Nx), density[Ny//4, :])
                axes[0].set_xlabel('X Position')
                axes[0].set_ylabel(f'Density ρ (x,y = {Ny//4})')
                save_path = os.path.join(path, f'density_decay_{step}.png')
                axes[0].grid()  
                figs[0].savefig(save_path,  pad_inches=1)
                

            else:                
                # 2. Save the sin flow of velocity at a point 
                
                axes[0].set_ylim([-eps, eps])
                axes[0].plot(np.arange(Ny), velocity[:, Nx//4, 1])
                axes[0].set_xlabel('Y Position')                    
                axes[0].set_ylabel(f'Velocity u(x = {Nx//4},y)')
                axes[0].grid()  
                
        
              
        
        """Calculate Theoretical Velocity"""
        if (exp_type == CV.velocity):
            kinematic_viscosity = 1/3*(1/omega-1/2)
            y_val = instant_theoretical_velocity(kinematic_viscosity)
            y_val = y_val*(v_not[:, Nx//4, 1])
            v_not[:, Nx//4, 1] = y_val
            theoretical_velocity.append(y_val.max())   

                
        """To get only the simulated maximas"""
        if exp_type == CV.velocity:
            den_or_vel_list.append(np.max(velocity[:, :, 1]))
            max_min_list.append(np.min(velocity[:, :, 1]))
            x_value.append(step)
        else:            
            den_or_vel.append(density[Ny//4,Nx//4]-1)
            den_or_vel_list.append(np.max(density - 1))
            max_min_list.append(np.min(density - 1))
            x_value.append(step)
        
    """View animation of density flow"""
    if(exp_type == CV.density):
        animate(density_animation)

    """View Sinusoidal Decay over total time"""
    if exp_type == CV.velocity:
        save_path = os.path.join(path, f'velocity_decay_{step}.png')
        figs[0].savefig(save_path,  pad_inches=1)   
        x = np.arange(steps)
        axes[2].cla()
        axes[2].set_xlim([0,len(x)])
            
        X_Y_Spline = make_interp_spline(np.arange(steps), den_or_vel_list)
        X_Y_Spline_1 = make_interp_spline(np.arange(steps), max_min_list)
        X_ = np.linspace(0,x_value[-1], 500)
        Y_ = X_Y_Spline(X_)
        Y_1 = X_Y_Spline_1(X_)
        axes[2].plot(X_, Y_, color='blue') 
        axes[2].plot(X_, Y_1, color='blue') 
        axes[2].plot(x, theoretical_velocity, color = 'black',linestyle='dotted',  linewidth=3)   
        axes[2].fill_between(X_,Y_,Y_1, color = "blue", alpha=.2)           
        axes[2].set_xlabel('Time evolution')
        axes[2].set_ylabel(f'Velocity at (y = {Ny//4})')
        axes[2].legend(
            ['Simulated Maxima','Simulated Minima', 'Analytical ux(y=25)']) 
        axes[2].grid()       
        save_path = os.path.join(path, f'omega_{omega}.png')
        figs[2].savefig(save_path,  pad_inches=1)
        
    else:
        den_or_vel_list = np.array(den_or_vel_list)
        x = argrelextrema(den_or_vel_list, np.greater)[0]
        den_or_vel_list = den_or_vel_list[x]
        
        x_value_t = np.arange(len(den_or_vel))        
        axes[2].plot(x_value_t, den_or_vel, color='blue')       
        axes[2].set_xlabel(f'Timestep (ω = {omega})')
        axes[2].set_ylabel(f'Density ρ (x = {Nx//4}, y = {Ny//4}')       
        
        save_path = os.path.join(path, f'omega_{omega}.png')
        axes[2].grid()  
        figs[2].savefig(save_path, bbox_inches='tight', pad_inches=0.5)

    """Kinematic Viscosity Calculation Analytical vs Simulation via curve fit"""
    
    simulated_viscosity = curve_fit(decay_perturbation, xdata=x, ydata=den_or_vel_list)[0][0]
    analytical_viscosity = (1/3) * ((1/omega) - 0.5)

    return simulated_viscosity, analytical_viscosity
