from matplotlib.lines import Line2D
import numpy as np
from lbm_common import lbm 

def couette_bounce_back(f,lid_vel,velocity):    

    """Bounce back and lid velocity exerted here"""
    rho_wall = (2 * (f[-1, :, 6] + f[-1, :, 2] + f[-1, :, 5]) + f[-1, :, 3] + f[-1, :, 0] + f[-1, :, 1])/(1+velocity[-1,:,1])
    
    # for bottom y = 0
    f[1, :, 2] = f[0, :, 4]
    f[1, :, 5] = f[0, :, 7]
    f[1, :, 6] = f[0, :, 8]
    # f[0, :, 4] = 0
    # f[0, :, 7] = 0
    # f[0, :, 8] = 0
    
    # for top y = max_size_y
    f[-2, :, 4] = f[-1, :, 2] 
    f[-2, :, 7] = f[-1, :, 5] - 1/6 *    lid_vel 
    # f[-2, :, 7] = f[-1, :, 5] - 1/2 *  rho_wall *  velocity[-1,:,1] 
    f[-2, :, 8] = f[-1, :, 6] + 1/6 *  lid_vel        
    # f[-2, :, 8] = f[-1, :, 6] + 1/2 *  rho_wall *  velocity[-1,:,1]        
    # f[-1, :, 2] = 0
    # f[-1, :, 5] = 0
    # f[-1, :, 6] = 0

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