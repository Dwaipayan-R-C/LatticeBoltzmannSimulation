"""Imports"""
from matplotlib.lines import Line2D
import numpy as np
from lbm_common import lbm 


def f_moving_wall(f, lid_vel):
    rho_wall = (2 * (f[-1, 1:-1, 6] + f[-1, 1:-1, 2] + f[-1, 1:-1, 5]) + f[-1, 1:-1, 3] + f[-1, 1:-1, 0] + f[-1, 1:-1, 1])    
    f[-2,1:-1,4] = f[-1,1:-1,2]
    f[-2,1:-1,7] = f[-1,1:-1,5] - 1/6  * rho_wall * lid_vel
    f[-2,1:-1,8] = f[-1,1:-1,6] + 1/6  * rho_wall*  lid_vel

    return f
    
def f_rigid_wall(f, top, down, left, right):

    if top:
        f[-2,:, [7,4,8]] = f[-1,:, [5,2,6]]
    if down:
        f[1,:, [5,2,6]] = f[0,:, [7,4,8]]
    if left:
        f[:,1, [8,1,5]] = f[:,0, [6,3,7]]
    if right:
        f[:,-2, [7,3,6]] = f[:,-1, [5,1,8]]
    return f
