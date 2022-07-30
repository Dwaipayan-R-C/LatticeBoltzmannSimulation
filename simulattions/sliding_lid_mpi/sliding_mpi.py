"""Library Imports"""
import os
import sys
import glob
import numpy as np
from mpi4py import MPI

#opposite directions for bouncing back
opposite_direction = np.array([0,3,4,1,2,7,8,5,6], dtype = np.int8)

# velocity constants for each particle 
W = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36])
c_x = np.array([ 0, 1, 0,-1, 0, 1,-1,-1, 1])
c_y = np.array([ 0, 0, 1, 0,-1, 1, 1,-1,-1])
c = np.array(list(zip(c_y,c_x)))


def density_calculation(f):
    """Calculates the density"""    
    density = np.sum(f, axis=2)
    return density


def periodic_boundary_with_pressure_variations(grid,rho_in,rho_out):
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
    velocity = (np.dot(f, c).T / density.T).T    
    return velocity

def streaming(f):
    """Streaming takes place here"""    
    for i in range(9):
        f[:,:,i] = np.roll(f[:,:,i],c[i], axis = (0,1))         
    return f


def calculate_collision(f, relaxation):
    """Collision calculated here"""   
    density = np.sum(f, axis=-1)    
    velocity = calculate_velocity(f, density)
    f_eq = calculate_equilibrium(density,velocity)
    f -= relaxation * (f-f_eq)    
    return f, density, velocity

def save_mpiio(comm, fn, g_kl):
    """
    Write a global two-dimensional array to a single file in the npy format
    using MPI I/O: https://docs.scipy.org/doc/numpy/neps/npy-format.html

    Arrays written with this function can be read with numpy.load.

    Parameters
    ----------
    comm
        MPI communicator.
    fn : str
        File name.
    g_kl : array_like
        Portion of the array on this MPI processes. This needs to be a
        two-dimensional array.
    """
    from numpy.lib.format import dtype_to_descr, magic
    magic_str = magic(1, 0)

    local_nx, local_ny = g_kl.shape
    nx = np.empty_like(local_nx)
    ny = np.empty_like(local_ny)

    commx = comm.Sub((True, False))
    commy = comm.Sub((False, True))
    commx.Allreduce(np.asarray(local_nx), nx)
    commy.Allreduce(np.asarray(local_ny), ny)

    arr_dict_str = str({ 'descr': dtype_to_descr(g_kl.dtype),
                         'fortran_order': False,
                         'shape': (np.asscalar(nx), np.asscalar(ny)) })
    while (len(arr_dict_str) + len(magic_str) + 2) % 16 != 15:
        arr_dict_str += ' '
    arr_dict_str += '\n'
    header_len = len(arr_dict_str) + len(magic_str) + 2

    offsetx = np.zeros_like(local_nx)
    commx.Exscan(np.asarray(ny*local_nx), offsetx)
    offsety = np.zeros_like(local_ny)
    commy.Exscan(np.asarray(local_ny), offsety)

    file = MPI.File.Open(comm, fn, MPI.MODE_CREATE | MPI.MODE_WRONLY)
    if comm.Get_rank() == 0:
        file.Write(magic_str)
        file.Write(np.int16(len(arr_dict_str)))
        file.Write(arr_dict_str.encode('latin-1'))
    mpitype = MPI._typedict[g_kl.dtype.char]
    filetype = mpitype.Create_vector(g_kl.shape[0], g_kl.shape[1], ny)
    filetype.Commit()
    file.Set_view(header_len + (offsety+offsetx)*mpitype.Get_size(),
                  filetype=filetype)
    file.Write_all(g_kl.copy())
    filetype.Free()
    file.Close()


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

def calculate_equilibrium(density, velocity, simulation = False):
    """Calculates the collision equlibrium function Feq"""
    if(simulation == True):
        vel_x2_y2 = velocity[:,0] ** 2 + velocity[:,1] ** 2
        cu = np.dot(velocity,c.T)
        squared_velocity = cu ** 2
        f_eq = ((1 + 3*(cu.T) + 9/2*(squared_velocity.T) - 3/2*(vel_x2_y2.T)) * density ).T * W
    else:
        vel_x2_y2 = velocity[:,:,0] ** 2 + velocity[:,:,1] ** 2
        cu = np.dot(velocity,c.T)
        squared_velocity = cu ** 2
        f_eq = ((1 + 3*(cu.T) + 9/2*(squared_velocity.T) - 3/2*(vel_x2_y2.T)) * density.T ).T * W
    
    return f_eq


def f_comm1(f, commCart):
   
    top_src, top_dst = commCart.Shift(0, -1)
    bot_src, bot_dst = commCart.Shift(0, +1)
    lef_src, lef_dst = commCart.Shift(1, -1)
    rig_src, rig_dst = commCart.Shift(1, +1)

    # for 4 in
    '''
    6 7 8
    3 4 5
    0 1 2 <- moving lid here
    '''
    
    p1 = f[  1,  :,:].copy()
    p2 = f[ -1,  :,:].copy()
    commCart.Sendrecv(p1, top_dst, recvbuf=p2, source=top_src) # dest 7  src 1
    f[ -1,  :,:] = p2
    
    p1 = f[ -2,  :,:].copy()
    p2 = f[  0,  :,:].copy()
    commCart.Sendrecv(p1, bot_dst, recvbuf=p2, source=bot_src) # dest 1 src 7
    f[  0,  :,:] = p2
    
    p1 = f[  :,  1,:].copy()
    p2 = f[  :, -1,:].copy()
    commCart.Sendrecv(p1, lef_dst, recvbuf=p2, source=lef_src) # dest 3 src 5
    f[  :, -1,:] = p2
    
    p1 = f[  :, -2,:].copy()
    p2 = f[ :,  0, :,].copy()
    commCart.Sendrecv(p1, rig_dst, recvbuf=p2, source=rig_src) # dst 5 src 3
    f[  :,  0, :] = p2
    
    return f

def f_wall_parallel(f, coords, i):
   
    if (coords[0] == 0):
        f = f_rigid_wall(f, False, True, False, False)
        # if(i == 10):
        #     print(f"the rank for coord00 is {rank}")
        
    if (coords[0] == (y_decomp - 1)):
        f = f_moving_wall(f, lid_vel)
        # if(i == 10):
        #     print(f"the rank for coord0y is {rank}")
        
    if (coords[1] == 0):
        f = f_rigid_wall(f, False, False, True, False)
        # if(i == 10):
        #     print(f"the rank for coord11 is {rank}")
        
    if (coords[1] == (x_decomp - 1)):
        f = f_rigid_wall(f, False, False, False, True)
        # if(i == 10):
        #     print(f"the rank for coord1x is {rank}")    
    return f

def sliding_lid_simulation(Nx: int, Ny: int, omega: float, save_every, steps, lid_vel, y_decomp, x_decomp):
    """ Calculates the sliding_lid flow for Nx by Ny D2Q9 lattice"""

    # omega = (2*re)/(6*Nx*lid_vel+re)          
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    CommCart = comm.Create_cart((y_decomp, x_decomp), periods=(False, False), reorder=False)    
    
    print("The rank is: {} and the coordinate is: {} ".format(rank, CommCart.Get_coords(rank)))
    coords = CommCart.Get_coords(rank)
    # """ Calculates the sliding_lid flow for Nx by Ny D2Q9 lattice"""    

    # # """Starts from here"""
    # # # initilization of the grids used    
    rows = Ny // y_decomp + 2
    columns = Nx // x_decomp + 2
    density = np.ones((rows,columns))
    velocity = np.zeros((rows,columns,2))
    f = calculate_equilibrium(density,velocity)
    
        
    # # Iteration starts from here
    for step in range(steps):
        print(f'{step+1}//{steps}', end="\r")  

        # Streaming, Bounceback and Collision
        f = f_comm1(f, CommCart)
        f = streaming(f)
        f = f_wall_parallel(f, coords, step)        
        f, density, velocity = calculate_collision(f, omega)    
        
        if (step%save_every == 0):
            save_mpiio(CommCart, f'data/ux_{step}.npy', velocity[:,:,1])
            save_mpiio(CommCart, f'data/uy_{step}.npy', velocity[:,:,0])

      
        
steps = int(sys.argv[1])#10000
save_every = int(sys.argv[2])#1000
length = int(sys.argv[3])#100
width = int(sys.argv[4])#100
omega = float(sys.argv[5])#1000
lid_vel = float(sys.argv[6])#0.1
y_decomp = int(sys.argv[7])#1
x_decomp = int(sys.argv[8])#1

files = glob.glob('data/*')
for j in files:
    try:
        os.remove(j)
    except:
        pass        

sliding_lid_simulation(length, width, omega, save_every, steps, lid_vel, y_decomp, x_decomp)
    
