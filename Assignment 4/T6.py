### Adsorbed surface relaxation
from ase.optimize import BFGS
from gpaw import GPAW
from ase.io import read
import numpy as np

metals     = ['Au','Pt','Rh']
adsorbants = ['O', 'CO']

calcs = ['Data/T3/Au.gpw',
         'Data/T3/Pt.gpw',
         'Data/T3/Rh.gpw'] # Calculator wavefunctions from T3

E = np.zeros((3,2))

for i in range(3):
    for j in range(2):
        # Surfaces with adsorbates
        surface = read(f'Data/Adsorbed/{metals[i]}_{adsorbants[j]}.xyz')

        calc = GPAW(calcs[i]) # Restart GPAW calcs from T3
        surface.set_calculator(calc)
        dyn = BFGS(surface,
                   trajectory=f'{metals[i]}_{adsorbants[j]}.traj',
                   logfile=f'{metals[i]}_{adsorbants[j]}.log')
        dyn.run(fmax=.01)

        E[i,j] = surface.get_potential_energy()
    
np.savetxt('adsorbed potential energies.txt',E)

### Energy calculations
# The orders are [Au,Pt,Rh] and [CO,O2]
E_surface = np.loadtxt('Data/T3/potential energies.txt') # Surface w/o adsorbate
E_gases   = np.loadtxt('Data/T4/gas energies.txt')
E_O,E_CO  = np.loadtxt('Data/T6/adsorbed potential energies.txt',
                        unpack=True) # Surfaces w adsorbate

# Finding adsorption energies (E_ad = E_constituents - E_tot)
E_ad_O  = E_surface + E_gases[0]/2 - E_O # E_gases is for O2, we approximate half for O
E_ad_CO = E_surface + E_gases[1] - E_CO 

# Activation energies
E_activation = .22 - .3*(E_ad_O + E_ad_CO)

print('Adsorption Energies ([Au,Pt,Rh])')
print(f'O:  {E_ad_O} eV')
print(f'CO: {E_ad_CO} eV','\n')

print('Activation energies ([Au,Pt,Rh])')
print(f'{E_activation} eV')