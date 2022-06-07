import numpy as np
from lbm_common import lbm 

def couette_bounce_back(f,lid_vel,velocity):    

    """Bounce back and lid velocity exerted here"""
    rho_wall = (2 * (f[-1, :, 6] + f[-1, :, 2] + f[-1, :, 5]) + f[-1, :, 3] + f[-1, :, 0] + f[-1, :, 1])/(1+velocity[-1,:,1])
    max_size_x = f.shape[1]-1  # x
    max_size_y = f.shape[0]-1  # y
    velocity = np.zeros((f.shape[0],f.shape[1],2))
    velocity[-1,:,1] = lid_vel
    # for bottom y = 0
    f[1, :, 2] = f[0, :, 4]
    f[1, :, 5] = f[0, :, 7]
    f[1, :, 6] = f[0, :, 8]
    f[0, :, 4] = 0
    f[0, :, 7] = 0
    f[0, :, 8] = 0
    
    # for top y = max_size_y
    f[max_size_y - 1, :, 4] = f[max_size_y, :, 2] 
    f[max_size_y - 1, :, 7] = f[max_size_y, :, 5] - 1/2 *  rho_wall *  velocity[-1,:,1] #2/5
    f[max_size_y - 1, :, 8] = f[max_size_y, :, 6] + 1/2 *  rho_wall *  velocity[-1,:,1]        
    f[max_size_y, :, 2] = 0
    f[max_size_y, :, 5] = 0
    f[max_size_y, :, 6] = 0

    return f


def poiseuille_bounce_back(f,lid_vel):
    # baunce back without any velocity gain
    # btw why am i doing this from 1:-1 ???
    # rho_wall = np.average(np.sum(f,axis=0)) # doesnt do anything
    # rho_wall = f[0,:,-2] + f[1,:,-2] + 2* f[2,:,2]+ f[3,:,-2] +2*f[5,:,-2]+2*f[6,:,-2] # neither do you
    # for bottom y = 0
    f[1, 1:-1, 2] = f[0, 1:-1, 4]
    f[1, 1:-1, 5] = f[0, 1:-1, 7]
    f[1, 1:-1, 6] = f[0, 1:-1, 8]
    # for top y = -1
    f[-2, 1:-1, 4] = f[-1, 1:-1, 2]
    # f[7, 1:-1, -2] = f[5, 1:-1, -1] - 1 / 6 * lid_vel * rho_wall[1:-1]
    # f[8, 1:-1, -2] = f[6, 1:-1, -1] + 1 / 6 * lid_vel * rho_wall[1:-1]
    # code doesnt do anything so idk
    f[-2, 1:-1, 7] = f[-1, 1:-1, 5] - 1 / 6 * lid_vel
    f[-2, 1:-1, 8] = f[-1, 1:-1, 6] + 1 / 6 * lid_vel

    return f

def periodic_boundary_with_pressure_variations(f,rho_in,rho_out):
    # get all the values
    # p_diff = 0.5*(rho_in - rho_out)
    density = lbm.density_calculation(f)
    velocity = lbm.calculate_velocity(f, density)
    # velocity[:,:,1] = velocity[:,:,1] + p_diff/
    equilibrium = lbm.calculate_equilibrium(density, velocity)    
    equilibrium_in = lbm.calculate_equilibrium_point(rho_in, velocity[:,-2,:])    
    f[:, 0, :] = equilibrium_in + (f[:, -2,:] - equilibrium[:, -2,:])
    equilibrium_out = lbm.calculate_equilibrium_point(rho_out, velocity[:,1,:])
    f[:, -1, :] = equilibrium_out + (f[:, 1, :] - equilibrium[:, 1, :])
    return f, density, velocity