from ase.io.trajectory import Trajectory
from ase.io import read, write
traj = Trajectory('nequipOpt.traj')
atoms = traj[-1]
write('POSCAR',atoms)
