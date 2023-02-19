from ase.optimize import BFGS
from ase.build import fcc111
from gpaw import GPAW, PW
from ase.optimize import GPMin
import numpy as np

lattice_constants = [4.179, 3.969, 3.839]
materialList      = ["Au","Pt", "Rh"]

k = (4,4,1)
cutoff = 450

E = np.zeros(3)

for i in range(3):
  material = materialList[i]
  surface = fcc111(material, (3, 3, 3), a=lattice_constants[i], vacuum=6.0)
  calc = GPAW(xc='PBE',
            mode=PW(cutoff),
            kpts=k)
  surface.set_calculator(calc)
  dyn = GPMin(surface)
  dyn.run(fmax=0.01, steps=100)
  surface.set_calculator(calc)

  E[i] = surface.get_potential_energy()

np.savetxt('energies.txt',E)

print('Surfaces energies:')

for i in range(3):
    print(f'{materialList[i]} = {E[i]} eV')