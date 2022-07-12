"""Imports"""
from matplotlib.lines import Line2D
import numpy as np
from lbm_common import lbm 

def couette_bounce_back(f,lid_vel,velocity):    

    """ Couette bounce back and lid velocity exerted here"""   
    # for bottom wall
    f[1, :, 2] = f[0, :, 4]
    f[1, :, 5] = f[0, :, 7]
    f[1, :, 6] = f[0, :, 8]
    
    # for top wall
    f[-2, :, 4] = f[-1, :, 2] 
    f[-2, :, 7] = f[-1, :, 5] - 1/6 *    lid_vel
    f[-2, :, 8] = f[-1, :, 6] + 1/6 *  lid_vel
    return f


def poiseuille_bounce_back(f,lid_vel): 
    """ Poiseuille bounce back exerted here"""   
    # Bottom wall
    f[1, 1:-1, 2] = f[0, 1:-1, 4]
    f[1, 1:-1, 5] = f[0, 1:-1, 7]
    f[1, 1:-1, 6] = f[0, 1:-1, 8]

    # Top wall
    f[-2, 1:-1, 4] = f[-1, 1:-1, 2]    
    f[-2, 1:-1, 7] = f[-1, 1:-1, 5] 
    f[-2, 1:-1, 8] = f[-1, 1:-1, 6] 
    return f



def sliding_bounce_back(f,lid_vel,velocity):    

    """Bounce back and lid velocity exerted here"""
    rho_wall = (2 * (f[-1, 1:-1, 6] + f[-1, 1:-1, 2] + f[-1, 1:-1, 5]) + f[-1, 1:-1, 3] + f[-1, 1:-1, 0] + f[-1, 1:-1, 1])/(1+velocity[-1,1:-1,1])
   
    # Bottom wall
    f[1,1:-1,2] = f[0,1:-1,4]
    f[1,1:-1,5] = f[0,1:-1,7]
    f[1,1:-1,6] = f[0,1:-1,8]

    # Top wall
    f[-2,1:-1,4] = f[-1,1:-1,2]
    f[-2,1:-1,7] = f[-1,1:-1,5] - 1/2 *  rho_wall *  lid_vel
    f[-2,1:-1,8] = f[-1,1:-1,6] + 1/2 *  rho_wall *  lid_vel

    # Right wall
    f[1:-1,1,1] = f[1:-1,0,3]
    f[1:-1,1,5] = f[1:-1,0,7]
    f[1:-1,1,8] = f[1:-1,0,6]

    # left wall
    f[1:-1,-2,3] = f[1:-1,-1,1]
    f[1:-1,-2,6] = f[1:-1,-1,8]
    f[1:-1,-2,7] = f[1:-1,-1,5]
    return f