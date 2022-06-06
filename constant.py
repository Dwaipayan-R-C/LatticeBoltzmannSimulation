import numpy as np

#opposite directions for bouncing back
opposite_direction = np.array([0,3,4,1,2,7,8,5,6], dtype = np.int8)

# weight values for each direction
W = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36])
            #    0, 1, 2, 3, 4, 5, 6, 7, 8
c_x = np.array([ 0, 1, 0,-1, 0, 1,-1,-1, 1])
c_y = np.array([ 0, 0, 1, 0,-1, 1, 1,-1,-1])
c = np.array(list(zip(c_y,c_x)))


"""Shear Wave constants"""
velocity = 'velocity'
density = 'density'

"""Couette Flow"""
couette_flow = 'couette_flow'

"""Experiment Global settings"""
shear_wave = 'shear_wave'

