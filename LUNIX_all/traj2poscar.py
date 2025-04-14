from ase.io import read, write

# 读取traj文件
atoms = read("C:/Users/renyq/Desktop/result/nequipOpt.traj")

# 写入POSCAR文件
write('POSCAR', atoms, format='vasp')