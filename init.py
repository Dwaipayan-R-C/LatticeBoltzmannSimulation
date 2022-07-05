"""Imports"""
from lbm_common import constant as CV
from simulattions import shear_wave as shear_wave
from simulattions import couette_flow as couette
from simulattions import poiseuilleFlow as poiseuille
from simulattions import sliding_lid as sliding_lid


def experiment_type(type):
    if(type == CV.shear_wave):
        """Global settings"""
        Nx = 50
        Ny = 50
        omega = 1.5
        eps = 0.1
        output_dir = 'results'
        save_every = 500
        simulation_type = 'velocity' # density/ velocity
        steps = 7000        
        simulated, analytical = shear_wave.shear_wave_simulation(Nx,Ny,omega,eps,output_dir,save_every,steps,simulation_type)
        print('Analytical Viscosity : '+ str(analytical))
        print('Simulated Viscosity : '+ str(simulated))
    
    elif (type==CV.couette_flow):        
        Nx = 100
        Ny = 50
        steps = 20000
        omega = 1.5
        output_dir = 'results'
        save_every = 500
        lid_velocity = 0.01
        couette.couette_flow_simulation(Nx,Ny,omega,output_dir,save_every,steps, lid_velocity)

    elif (type == CV.poiseuille_flow):
        Nx = 50
        Ny = 50
        steps = 10000
        rho_null = 1
        p_diff = 0.001
        omega = 1
        output_dir = 'results'
        save_every = 200        
        poiseuille.poiseuille_simulation(rho_null,p_diff, output_dir, Nx,Ny,omega,steps, save_every)

    elif (type == CV.sliding_lid):
        Nx = 300
        Ny = 300
        steps = 20001
        rho_null = 1
        re = 1000
        output_dir = 'results'
        save_every = 500    
        lid_velocity = 0.1    
        sliding_lid.sliding_lid_simulation(Nx, Ny, re, output_dir, save_every, steps, lid_velocity)

"""Define the experiemnt  category"""
if __name__ == "__main__":
    # 'shear_wave', 'couette_flow', 'poiseuille_flow', 'sliding_lid'
    experiment_type("sliding_lid")

# 