from ase.build import molecule
from gpaw import PW,GPAW
import numpy as np

moleculeList = ["CO", "O2"] 
cutoff = 450 # eV
roomtemp = 300 # k
atmosphere = 10**5 #bar

E = np.zeros(2)

for i in range(2):
    atoms = molecule(moleculeList[i], cell=(12, 12, 12))
    atoms.center()
    calc = GPAW(xc='PBE',
                mode=PW(cutoff),
                kpts={'gamma': True})
    atoms.set_calculator(calc)

    E[i] = atoms.get_potential_energy()

np.savetxt('Data/T4/gas energies.txt',E)