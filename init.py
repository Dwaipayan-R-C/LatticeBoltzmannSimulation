"""Imports"""

from lbm_common import constant as CV
from simulattions import shear_wave as shear_wave
from simulattions import couette_flow as couette
from simulattions import poiseuilleFlow as poiseuille
from simulattions import sliding_lid as sliding_lid
import os
import matplotlib.pyplot as plt
import numpy as np
# import ipyparallel as ipp

# Experimental setup


def experiment_type(type):
    plt.rcParams.update({'font.size': CV.fontsize})
    plt.rcParams["font.family"] = CV.fontfamily

    # Shear wave settings
    if(type == CV.shear_wave):
        """Global settings"""
        # These are the omega values taken to measure for Analytical vs Simulated kinematic viscosity
        # Turn on the kinematic_vs_analytical to True in order to get a result of the same.
        # In that case, please comment other plots in shear wave decay
        omega_list = [0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2]
        analytical_viscosity, simulated_viscosity = [], []
        Nx = 50
        Ny = 50
        omega = 1
        eps = 0.01
        output_dir = 'results'
        save_every = 100
        simulation_type = 'velocity'  # density/ velocity
        steps = 2000
        # Please turn off all the plots inside shear wave decay
        kinematic_vs_analytical = False

        # region kinematic viscosity vs analytical viscosity
        if(kinematic_vs_analytical == True and simulation_type == "velocity"):
            for i in omega_list:
                simulated, analytical = shear_wave.shear_wave_simulation(
                    Nx, Ny, i, eps, output_dir, save_every, steps, simulation_type)
                analytical_viscosity.append(analytical)
                simulated_viscosity.append(simulated)

            common_path = os.path.join(output_dir, 'shear_decay')
            fig, ax = plt.subplots()
            ax.plot(np.array(omega_list), np.array(simulated_viscosity))
            ax.plot(np.array(omega_list), np.array(analytical_viscosity),
                    color='r', linestyle='dotted',  linewidth=3)
            ax.grid()
            ax.set_xlabel("Omega")
            ax.set_ylabel("Viscosity")
            ax.set_title("Kinematic viscosity vs Omega")
            ax.legend(['Simulated viscosity', 'Analytical viscosity'])
            save_path = os.path.join(common_path, f'viscosity_vs_omega.png')
            fig.savefig(save_path, pad_inches=1)
        # endregion
        else:
            simulated, analytical = shear_wave.shear_wave_simulation(
                Nx, Ny, omega, eps, output_dir, save_every, steps, simulation_type)
            print('Analytical Viscosity : ' + str(analytical))
            print('Simulated Viscosity : ' + str(simulated))
    # Couette flow settings
    elif (type == CV.couette_flow):
        Nx = 100
        Ny = 50
        steps = 20000
        omega = 1.5
        output_dir = 'results'
        save_every = 1000
        lid_velocity = 0.1
        couette.couette_flow_simulation(
            Nx, Ny, omega, output_dir, save_every, steps, lid_velocity)

    # Poisuille flow settings
    elif (type == CV.poiseuille_flow):
        Nx = 50
        Ny = 50
        steps = 10000
        rho_null = 1
        p_diff = 0.001
        omega = 1
        output_dir = 'results'
        save_every = 500
        poiseuille.poiseuille_simulation(
            rho_null, p_diff, output_dir, Nx, Ny, omega, steps, save_every)

    # Sliding lid serialized settings
    elif (type == CV.sliding_lid):
        Nx = 150
        Ny = 150
        steps = 10001
        rho_null = 1
        re = 300
        output_dir = 'results'
        save_every = 2000
        lid_velocity = 0.1        
        sliding_lid.sliding_lid_simulation(Nx, Ny, re, output_dir, save_every, steps, lid_velocity)

"""Define the experiemnt  category"""
if __name__ == "__main__":
    # 'shear_wave', 'couette_flow', 'poiseuille_flow', 'sliding_lid'
    experiment_type("sliding_lid")

#
