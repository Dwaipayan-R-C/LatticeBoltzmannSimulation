import numpy as np
import constant as CV
import shear_wave as shear_wave
import couette_flow as couette

def experiment_type(type):
    if(type == CV.shear_wave):
        """Global settings"""
        Nx = 150
        Ny = 100
        omega = 1
        eps = 0.01
        output_dir = 'results'
        save_every = 100
        simulation_type = 'density' # density/ velocity
        steps = 1000
        simulated, analytical = shear_wave.shear_wave_simulation(Nx,Ny,omega,eps,output_dir,save_every,steps,simulation_type)
        print('Analytical Viscosity : '+ str(analytical))
        print('Simulated Viscosity : '+ str(simulated))
    
    elif (type==CV.couette_flow):        
        Nx = 100
        Ny = 50
        steps = 2000
        omega = 1
        output_dir = 'results'
        save_every = 200
        lid_velocity = 1
        couette.couette_flow_simulation(Nx,Ny,omega,output_dir,save_every,steps, lid_velocity)

"""Define the experiemnt  category"""
if __name__ == "__main__":
    # 'shear_wave', 'couette_flow'
    experiment_type("shear_wave")
