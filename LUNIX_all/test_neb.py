from ase.optimize import BFGS, FIRE
from ase.io import read, write
from nequip.ase import NequIPCalculator
from ase.neb import NEB

model_path = '/public/home/ac877eihwp/haoyl/potential/test/again/prototypeModel.pth'
calc = NequIPCalculator.from_deployed_model(model_path, device='cpu')


def var_name(var, all_var=locals()):
    return [var_name for var_name in all_var if all_var[var_name] is var][0]


def run_neb(IS, FS, nImg, out_file):
    steps = 2000
    a = False
    while not a:
        a = cycle_neb(IS, FS, nImg, steps)
        nImg += 1

    write('%s.xyz' % out_file, a[:])

    TS = sorted(a, key=lambda x: x.get_potential_energy())[-1]

    return TS


def cycle_neb(IS, FS, nImg, steps):
    a = [IS.copy() for i in range(nImg+1)] + [FS]
    for i in a: 
        i.calc = NequIPCalculator.from_deployed_model(model_path, device='cpu')
    neb = NEB(a, climb=True)
    neb.interpolate(method='idpp', mic=True)
    if FIRE(neb).run(fmax=0.05, steps=steps):
        return a
    else:
        return False


def read_data(file_name):
    Atoms = read(file_name)
    Atoms.calc = calc
    BFGS(Atoms).run(fmax=0.01)
    return Atoms

IS, FS = read_data('POSCAR'), read_data('POSCAR2')
#write('IS7.xyz', IS)
#write('FS7.xyz', FS)

TS = run_neb(IS, FS, 8, 'TS')

d = [IS, TS, FS]

print('-'*50)
for i in d:
    i.calc = NequIPCalculator.from_deployed_model(model_path, device='cpu')
    print('%8s : %.3f' % (var_name(i), i.get_potential_energy()))
    
write('R7.xyz', d[:])
