from ase.optimize import BFGS
from ase.build import fcc111
from gpaw import GPAW, PW
from ase.optimize import GPMin, BFGS
from ase.io.trajectory import Trajectory
import numpy as np

lattice_constants = [4.179, 3.969, 3.839]
E_bulk            = [-3.146,-6.431,-7.306]
materialList      = ['Au','Pt','Rh']

k = (4,4,1)
cutoff = 450

E_pot  = np.zeros(3)
E_surf = np.zeros(3)
area   = np.zeros(3)

for i in range(3):
  material = materialList[i]
  surface = fcc111(material, (3, 3, 3), a=lattice_constants[i], vacuum=6.0)
  calc = GPAW(xc='PBE',
            mode=PW(cutoff),
            kpts=k)
  surface.set_calculator(calc)
  dyn = BFGS(surface,trajectory=f'{materialList[i]}.traj',
             logfile=f'{materialList[i]}.log')
  dyn.run(fmax=0.01, steps=100)

  calc.write(f'{materialList[i]}.gpw') # Save calc for coming tasks

  E_pot[i] = surface.get_potential_energy()
  
  # Finding area of the surface (A_parallelogram = ABsin(theta_AB))

  dims_angles = surface.get_cell_lengths_and_angles()
  area[i] = dims_angles[0]*dims_angles[1]*np.sin(dims_angles[5]*np.pi/180)

  # Calculating the surface energy acc to E_S = (E_slab - N*E_bulk)/2A
  E_surf[i] = (E_pot[i] - surface.get_number_of_atoms()*E_bulk[i]) / (2*area[i])

# Saving relevant energies
np.savetxt('Potential energies.txt',E_pot)
np.savetxt('Surface energies.txt',E_surf)
np.savetxt('Areas.txt',area)
print('Surfaces energies:')

for i in range(3):
    print(f'{materialList[i]} = {E_surf[i]} eV')