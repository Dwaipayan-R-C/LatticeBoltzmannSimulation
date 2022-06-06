import numpy as np
import constant as CV
import matplotlib.pyplot as plt

"""Basic lbm functions required to stream and collide"""
def density_calculation(f):
    """Calculates the density"""    
    density = np.sum(f, axis=2)
    return density

def calculate_velocity(f, density):
    """Calculates the velocity"""   
    velocity = np.dot(f, CV.c) / (density[:, :, None]) 
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
    return f, density,velocity

def bounce_back(f,lid_vel,velocity):    
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
    # f[max_size_y - 1, :, 7] = f[max_size_y, :, 5] - 0.5 * lid_vel  * rho_wall
    f[max_size_y - 1, :, 8] = f[max_size_y, :, 6] + 1/2 *  rho_wall *  velocity[-1,:,1]
    # f[max_size_y - 1, :, 8] = f[max_size_y, :, 6] + 0.5 * lid_vel  * rho_wall
        
    f[max_size_y, :, 2] = 0
    f[max_size_y, :, 5] = 0
    f[max_size_y, :, 6] = 0

    return f



def calculate_equilibrium(density, velocity):
    """Calculates the collision equlibrium function Feq"""
    cu = np.dot(velocity,CV.c.T)
    squared_velocity = cu ** 2
    vel_x2_y2 = velocity[:,:,0] ** 2 + velocity[:,:,1] ** 2
    f_eq = ((1 + 3*(cu.T) + 9/2*(squared_velocity.T) - 3/2*(vel_x2_y2.T)) * density.T ).T * CV.W
    return f_eq

def plot(velocity,ax=plt):
    """Calculates Velocity profile """
    v = np.sqrt(velocity[:, :, 0]**2 + velocity[:, :, 1]**2)
    ax.imshow(v, cmap='RdBu_r', vmin=0, interpolation='spline16')


