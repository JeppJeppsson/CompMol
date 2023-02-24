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
        surface = read(f'Data/Adsorbed/{metals[i]}_{adsorbants[j]}.xyz') # Surfaces with absorbates
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

# Calculate the adsorption energies according to
# E_adsorption = E_with_adsorbate - (E_without_adsorbate + E_gas/2)
E_O_ads  = E_O - (E_surface - E_gases[0]/4) # We approximate E_O = E_O2/2
E_CO_ads = E_CO - (E_surface - E_gases[1]/2)
O_adsorp  = E_surface-E_O
CO_adsorp = E_surface-E_CO

# Activation energies
E_activation = .22 - .3*(O_adsorp + CO_adsorp)

print('Adsorption Energies ([Au,Pt,Rh])')
print(f'O:  {O_adsorp} eV')
print(f'CO: {CO_adsorp} eV','\n')

print('Activation energies ([Au,Pt,Rh])')
print(f'{E_activation} eV')