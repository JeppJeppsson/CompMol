# Initialization
import numpy as np
import scipy as sp
import ase
from ase import Atoms
from ase.io.trajectory import Trajectory

# Importing the trajectory
traj = Trajectory('temp.traj')

print(traj[-1])

#for atoms in traj:
#