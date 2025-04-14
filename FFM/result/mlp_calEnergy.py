import ase.io
from nequip.ase import NequIPCalculator
from ase.optimize import BFGS
import ase.io
import os
import sys
import logging

nequipModel=NequIPCalculator.from_deployed_model(model_path='/public/home/ac877eihwp/renyq/val/prototypeModel.pth',device='cpu')

#vasp_energies = []
#nequip_energies = []
#ss=ase.io.read('OUTCAR', index=':')
#for i,s in enumerate(ss):
#  vaspE=s.get_potential_energy()
#  s.set_calculator(nequipModel)
#  nequipE=s.get_potential_energy()
#  log_message = f"{i}: vaspE = {vaspE}, nequipE = {nequipE}"
#  logging.info(log_message)
#  print(log_message)
# vasp_energies.append(vaspE)
#  nequip_energies.append(nequipE)

struct=ase.io.read('POSCAR')
struct.set_calculator(nequipModel)
print(f' Starting optmization by NequIP model:')
optJob=BFGS(struct, trajectory='nequipOpt.traj')
optJob.run(fmax=0.05,steps=1000)
