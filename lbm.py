import numpy as np
import constant as CV
import matplotlib.pyplot as plt

"""Basic lbm functions required to stream and collide"""
def density_calculation(f):
    """Calculates the density"""
    # should calculate whole grid
    density = np.sum(f, axis=2)
    return density

def calculate_velocity(f):
    """Calculates the velocity"""
    density = density_calculation(f)
    velocity = (np.dot(f, CV.c).T / density.T).T
    return velocity

def streaming(f):
    """Streaming takes place here"""
    for i in range(9):
        f[:, :, i] = np.roll(np.roll(f[:, :, i], CV.c[i, 0], axis=0), CV.c[i, 1], axis=1)    
    return f

# Calculate Feq
def calculate_equilibrium(density, velocity):
    """Calculates the collision equlibrium function Feq"""
    cu = np.dot(velocity,CV.c.T)
    squared_velocity = cu ** 2
    vel_x2_y2 = velocity[:,:,0] ** 2 + velocity[:,:,1] ** 2
    f_eq = ((1 + 3*(cu.T) + 9/2*(squared_velocity.T) - 3/2*(vel_x2_y2.T)) * density.T ).T * CV.W
    return f_eq

def collision(f, omega):
    """Calculates the collision"""
    # density
    density = density_calculation(f)
    # velocity
    u = calculate_velocity(f)
    # calculates feq
    feq = calculate_equilibrium(density, u)
    # calculates f after collision
    f += omega * (feq - f)
    return f, density, u

def plot(velocity,ax=plt):
    """Calculates Velocity profile """
    v = np.sqrt(velocity[:, :, 0]**2 + velocity[:, :, 1]**2)
    ax.imshow(v, cmap='RdBu_r', vmin=0, interpolation='spline16')


