# Initialization
import time
start_time = time.time()

import gpaw as gp
from gpaw import GPAW, PW
from ase import Atoms
from ase.io import read
from ase.optimize import BFGS
from ase.visualize import view

# Importing the atoms
deca = read('Na6-structures/christmas-tree.xyz')
tree = read('Na6-structures/half-decahedron.xyz')

# Setting up the GPAW calculator
calc_deca = gp.GPAW(mode=PW(200),xc='LDA',txt='pw200_LDA/pw200_LDA_deca.txt')
calc_tree = gp.GPAW(mode=PW(200),xc='LDA',txt='pw200_LDA/pw200_LDA_tree.txt')

deca.set_calculator(calc_deca)
tree.set_calculator(calc_tree)

# Optimizing
opt_deca = BFGS(deca,logfile='pw200_LDA/bfgs_deca.log')
opt_tree = BFGS(tree,logfile='pw200_LDA/bfgs_tree.log')

opt_deca.run(fmax=.05)
opt_tree.run(fmax=.05)

# Extracting the energy of the relaxed cluster
E_deca = deca.get_potential_energy()
E_tree = tree.get_potential_energy()

print(f'Ground state energy: {E_tree} eV')
print(f'Second to lowest state energy: {E_deca} eV')

print('Total time: {:1f} seconds'.format(time.time() - start_time))