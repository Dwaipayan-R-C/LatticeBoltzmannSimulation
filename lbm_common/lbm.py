"""Imports"""
import numpy as np
from lbm_common import constant as CV
import matplotlib.pyplot as plt
from lbm_common import boundary as boundary

"""Basic lbm functions required to stream and collide"""
def density_calculation(f):
    """Calculates the density"""    
    density = np.sum(f, axis=2)
    return density


def periodic_boundary_with_pressure_variations(grid,rho_in,rho_out):
    """For poiseuille PBC is required 

    Args:
        grid : F - The spacial grid
        rho_in : Inlet density
        rho_out : Outlet density

    Returns:
        density,velocity, f: density,velocity, f of the grid
    """
    rho = density_calculation(grid)
    velocity = calculate_velocity(grid,rho)
    # overall equilibrium
    equilibrium = calculate_equilibrium(rho,velocity)
    equilibrium_in = calculate_equilibrium(rho_in, velocity[:,-2,:], True)
    # equilibrium for inlet 1,5,8
    grid[:,0, :] = equilibrium_in + (grid[:,-2, :] - equilibrium[:,-2, :])
    # equilibrium for outlet 3,6,7
    equilibrium_out = calculate_equilibrium(rho_out, velocity[:, 1, :], True)    
    grid[:,-1, :] = equilibrium_out + (grid[:,1, :] - equilibrium[:,1, :])
    return rho, velocity, grid


def calculate_velocity(f, density):
    """Calculates the velocity"""   
    velocity = (np.dot(f, CV.c).T / density.T).T    
    return velocity

def streaming(f):
    """Streaming takes place here"""    
    for i in range(9):
        f[:,:,i] = np.roll(f[:,:,i],CV.c[i], axis = (0,1))         
    return f


def calculate_collision(f, relaxation):
    """Collision calculated here"""   
    density = np.sum(f, axis=-1)    
    velocity = calculate_velocity(f, density)
    f_eq = calculate_equilibrium(density,velocity)
    f -= relaxation * (f-f_eq)    
    return f, density, velocity


def calculate_equilibrium(density, velocity, simulation = False):
    """Calculates the collision equlibrium function Feq"""
    if(simulation == True):
        vel_x2_y2 = velocity[:,0] ** 2 + velocity[:,1] ** 2
        cu = np.dot(velocity,CV.c.T)
        squared_velocity = cu ** 2
        f_eq = ((1 + 3*(cu.T) + 9/2*(squared_velocity.T) - 3/2*(vel_x2_y2.T)) * density ).T * CV.W
    else:
        vel_x2_y2 = velocity[:,:,0] ** 2 + velocity[:,:,1] ** 2
        cu = np.dot(velocity,CV.c.T)
        squared_velocity = cu ** 2
        f_eq = ((1 + 3*(cu.T) + 9/2*(squared_velocity.T) - 3/2*(vel_x2_y2.T)) * density.T ).T * CV.W
    
    return f_eq

def plot(velocity,ax=plt):
    """Calculates Velocity profile """
    v = np.sqrt(velocity[:, :, 0]**2 + velocity[:, :, 1]**2)
    ax.imshow(v, cmap='RdBu_r', vmin=0, interpolation='spline16')


