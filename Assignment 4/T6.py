from ase.optimize import BFGS
from gpaw import GPAW
from ase.io import read
import numpy as np

metals     = ['Au','Pt','Rh']
adsorbants = ['O', 'CO']

calcs = ['Data/T3/Au.gpw',
         'Data/T3/Pt.gpw',
         'Data/T3/Rh.gpw']

E = np.zeros((3,2))

for i in range(3):
    for j in range(2):
        surface = read(f'Data/Adsorbed/{metals[i]}_{adsorbants[j]}.xyz') # Read ads. surfaces
        calc = GPAW(calcs[i]) # Restart GPAW cals from T3
        surface.set_calculator(calc)
        dyn = BFGS(surface,
                   trajectory=f'Data/T6/{metals[i]}_{adsorbants[j]}.traj',
                   logfile=f'Data/T6/{metals[i]}_{adsorbants[j]}.log')
        dyn.run(fmax=.01)
        E[i,j] = surface.get_potential_energy()
    
np.savetxt('energies Au Pt Rh O CO.txt',E)