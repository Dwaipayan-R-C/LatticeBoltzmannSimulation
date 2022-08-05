import numpy as np


#opposite directions for bouncing back
opposite_direction = np.array([0,3,4,1,2,7,8,5,6], dtype = np.int8)

# velocity constants for each particle 
W = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36])
c_x = np.array([ 0, 1, 0,-1, 0, 1,-1,-1, 1])
c_y = np.array([ 0, 0, 1, 0,-1, 1, 1,-1,-1])
c = np.array(list(zip(c_y,c_x)))


"""Experiment constants"""
velocity = 'velocity'
density = 'density'
couette_flow = 'couette_flow'
poiseuille_flow = 'poiseuille_flow'
shear_wave = 'shear_wave'
sliding_lid = 'sliding_lid'


"""Plot Global settings"""
fontsize = 11
fontfamily = 'Cambria'

