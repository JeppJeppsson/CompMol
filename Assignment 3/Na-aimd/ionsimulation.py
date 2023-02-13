from gpaw import GPAW
from ase import Atoms
from ase.io import read, Trajectory
from ase.units import fs,kB
from ase.md.npt import NPT

atoms = read('startconfig.xyz')

calc = GPAW(
    mode = 'lcao',
    xc = 'PBE',
    basis = 'dzp',
    symmetry = {'point_group':False},
    charge = 1,
    txt = 'output.gpaw-out')

atoms.set_calculator(calc)

dyn = NPT(
    atoms,
    temperature_K = 350,
    timestep = .5*fs, # This is not an appropriate timestep , I can tell you that !
    ttime = 20*fs, # Don ’t forget the fs!
    externalstress = 0, # We don ’t use the barostat , but this needs to be set anyway !
    logfile = 'mdOutput.log')

trajectory = Trajectory('someDynamics.traj','w',atoms)

dyn.attach(trajectory.write, interval=1) # Write the current positions etc. to file each timestep

dyn.run(4000) # Run 10 steps of MD simulation