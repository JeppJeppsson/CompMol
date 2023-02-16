# Initialization
import numpy as np
import scipy as sp
from ase import Atoms
from ase.io.trajectory import Trajectory
from matplotlib import pyplot as plt

# Importing the trajectory
traj = Trajectory('Dynamics.traj')

#Parameters
start_snap = 0
stop_snap = len(traj) 
distances = np.zeros([])
distance_sum = np.zeros(72)
cell_volume = traj[0].get_volume()
particle_density = cell_volume/len(traj[0]) 

for i in range(start_snap,stop_snap):
    atoms = traj[i]
    distances = np.append(distances,atoms.get_distances(72,atoms.get_atomic_numbers()==8,mic=True))
        
distance_sum *= (len(traj))**-1
distance_sum = np.sort(distance_sum)

plt.hist(distances,bins = 100)
hist = np.histogram(distances,bins=100)

dr = hist[1][1]-hist[1][0]
r = hist[1][1:]-dr/2
dn_r = hist[0]/(stop_snap-start_snap) 

gPrime = dn_r*(4*np.pi*r**2*dr*particle_density)**-1
g = dn_r*(dr)**-1

r_1 = 3
r_2 = 8
minapprox = int((r_1-r[0])/dr) #Place to start searching for first minimum
minimum = np.argmin(gPrime[minapprox-r_2:minapprox+r_2])+minapprox-r_2
plt.figure(figsize=(8, 6))
plt.plot(r[20:],gPrime[20:])

plt.xlabel("r [Å]",fontsize=20)
plt.ylabel("g(r)",fontsize=20)
plt.tight_layout()
plt.savefig('RadialDensityFunctionTask_2.pdf')
# print("solvation shell = ",np.trapz(g[:minimum],dx=dr))
print("solvation shell =  ", sum(hist[0][:minimum])/(stop_snap-start_snap))