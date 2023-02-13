from gpaw import GPAW
from ase.io import read, Trajectory
from ase.units import fs
from ase.md.npt import NPT

atoms = read('pepsi.xyz')

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
    timestep = .5*fs,
    ttime = 20*fs,
    externalstress = 0,
    logfile = 'mdOutput.log')

trajectory = Trajectory('someDynamics.traj','w',atoms)

dyn.attach(trajectory.write, interval=1)

dyn.run(4000)