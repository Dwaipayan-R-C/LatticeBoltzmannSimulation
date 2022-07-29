

import numpy as np
import matplotlib.pyplot as plt

lid_vel = 0.2
omega = 1.6
length = 100

calc_neu = lambda omega: (1/3 * (1/omega - 1/2)) 
calc_reynolds = lambda omega, length, vel: (vel * length) / calc_neu(omega)

ux = []
uy = []
                 
for i in range(0, 10000, 1000):
    var = 'sliding_lid_mpi/data/ux_{}.npy'.format(i)
    print(var)
    ux.append(np.load(var))
    uy.append(np.load('sliding_lid_mpi/data/uy_{}.npy'.format(i)))
    
ux = np.array(ux)
uy = np.array(uy)

print("the shape of acquired array is")
print(ux.shape)
nx, ny = ux[1,:,:].shape

i=1
fig=plt.figure(figsize=(20,8))
fig.suptitle("Sliding lid with omega = {} and lid vel = {} and reynold's number = {}".format(omega, lid_vel, calc_reynolds(omega, length, lid_vel)))

for itx in range(0, 10, 1):
        ax = plt.subplot(2,5,i)
        ax.set_title("iteration: {}".format(itx))
        ax.set_xlabel("length")
        ax.set_ylabel("width")
        ax.set_xlim(0,nx)
        ax.set_ylim(0,ny)
        #ax.invert_yaxis()
        intensity = ax.streamplot(np.arange(0, nx), np.arange(0,ny), ux[itx], uy[itx])
        i+=1

plt.tight_layout()
plt.savefig('plots/sliding_lid_parallelized.png')
plt.show()