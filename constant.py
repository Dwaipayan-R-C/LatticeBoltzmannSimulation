import numpy as np

#opposite directions for bouncing back
opposite_direction = np.array([0,3,4,1,2,7,8,5,6], dtype = np.int8)

# weight values for each direction
W = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36])
c_x = np.array([ 0, 1, 0,-1, 0, 1,-1,-1, 1])
c_y = np.array([ 0, 0, 1, 0,-1, 1, 1,-1,-1])
c = np.array(list(zip(c_x,c_y)))

"""Shear Wave constants"""
velocity = 'velocity'
density = 'density'

"""Experiment Global settings"""
shear_wave = 'shear_wave'

