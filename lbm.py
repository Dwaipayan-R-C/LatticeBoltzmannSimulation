import numpy as np
import constant as CV
import matplotlib.pyplot as plt

"""Basic lbm functions required to stream and collide"""
def density_calculation(f):
    # should calculate whole grid
    density = np.sum(f, axis=2)
    return density

def calculate_velocity(f):
    # first calculate rho for the given f
    density = density_calculation(f)
    velocity = (np.dot(f, CV.c).T / density.T).T
    return velocity

def streaming(f):
    # stream
    for i in range(9):
        f[:, :, i] = np.roll(np.roll(f[:, :, i], CV.c[i, 0], axis=0), CV.c[i, 1], axis=1)    
    return f

# Calculate Feq
def calculate_equilibrium(density, velocity):
    cu = np.dot(velocity,CV.c.T)
    squared_velocity = cu ** 2
    vel_x2_y2 = velocity[:,:,0] ** 2 + velocity[:,:,1] ** 2
    f_eq = ((1 + 3*(cu.T) + 9/2*(squared_velocity.T) - 3/2*(vel_x2_y2.T)) * density.T ).T * CV.W
    return f_eq

def collision(f, omega):
    density = density_calculation(f)
    u = calculate_velocity(f)
    # to compute the new f, first we need to compute feq
    feq = calculate_equilibrium(density, u)
    # relaxation
    f += omega * (feq - f)
    return f, density, u

# Calculates Velocity profile 
def plot(velocity,ax=plt):
    v = np.sqrt(velocity[:, :, 0]**2 + velocity[:, :, 1]**2)
    ax.imshow(v, cmap='RdBu_r', vmin=0, interpolation='spline16')


