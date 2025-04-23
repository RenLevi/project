'''v1.0.7'''
import ast
from ase import Atoms
from ase.build import rotate
from ase.build import hcp0001
from ase.io import read, write
import numpy as np
from ase.constraints import FixAtoms
import copy
import os
import re
#str to code
def char_to_value_extended(c):
    """将字符转换为数值（支持字母、数字和常见符号）"""
    if '0' <= c <= '9':
        return ord(c) - ord('0')  # 0-9
    elif 'a' <= c <= 'z':
        return ord(c) - ord('a') + 10  # 10-35
    elif 'A' <= c <= 'Z':
        return ord(c) - ord('A') + 36  # 36-61
    else:
        # 常见符号映射到62-93
        common_symbols = " !@#$%^&*()_+-=[]{};':\",./<>?\\|`~"
        if c in common_symbols:
            return 62 + common_symbols.index(c)
        else:
            raise ValueError(f"不支持的字符: {c}")

def value_to_char_extended(v):
    """将数值转换回字符"""
    if 0 <= v <= 9:
        return chr(v + ord('0'))
    elif 10 <= v <= 35:
        return chr(v - 10 + ord('a'))
    elif 36 <= v <= 61:
        return chr(v - 36 + ord('A'))
    else:
        common_symbols = " !@#$%^&*()_+-=[]{};':\",./<>?\\|`~"
        if 62 <= v < 62 + len(common_symbols):
            return common_symbols[v - 62]
        else:
            raise ValueError(f"无效值: {v}")

def encode_extended(s):
    """编码包含字母、数字和符号的字符串"""
    num = 0
    for c in s:
        num = num * 94 + char_to_value_extended(c)  # 总共94种字符
    return num

def decode_extended(num):
    """解码扩展编码的数字"""
    if num == 0:
        return "0"
    
    s = []
    while num > 0:
        s.append(value_to_char_extended(num % 94))
        num = num // 94
    return ''.join(reversed(s))


## 模型保存
def save_file(save_path,filename,model):
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    write(os.path.join(save_path, filename), model, format='vasp', vasp5=True, direct=True)
    with open(os.path.join(save_path, filename), 'r') as f:
        lines = f.readlines()
        print("".join(lines[5:8]))  # 查看坐标行示例
## 建模金属表面
def build_slab(element,size,vacuum,fix_num,save_path):#元素，尺寸，真空层厚度，保存格式，保存路径:mol_to_ad，固定底层原子层数
    x=size[0]
    y=size[1]
    #z=size[2]
    matel_surface = hcp0001(element, size=size, vacuum=vacuum)#可更改表面
    z_coords = matel_surface.positions[:, 2]
    threshold = fix_num*x*y
    fixed_indices = list(range(0, threshold))
    matel_surface.constraints = [FixAtoms(indices=fixed_indices)]
    slabname = element+'_hcp0001'+'.vasp'
    slab_save_path = save_path+'/slab'#保存至子文件夹
    save_file(slab_save_path,slabname,matel_surface)
    return matel_surface
## 分子操作
'''rotate'''
def rotate_mol(mol,rotate_matrix1,rotate_matrix2):#rotate molecule
    molcopy = copy.deepcopy(mol)
    (angle1,axis1) = rotate_matrix1
    (angle2,axis2) = rotate_matrix2
    molcopy.rotate(angle1,axis1,center=(0,0,0))
    molcopy.rotate(angle2,axis2,center=(0,0,0))
    return molcopy
'''place'''
def place_mol_on_surface(mol,surface,shift_vector):#place molecule
    surfacecopy = copy.deepcopy(surface)
    # 找到Ru(0001)表面最上层原子的z坐标最大值
    z_max = max(surfacecopy.positions[:, 2])
    # 计算分子的结合位点,将分子的质心移动到这个高度
    molecule_center = shift_vector + np.array([0,0,z_max])
    coordinates = mol.get_positions()
    average_coordinates = coordinates.mean(axis=0)
    mol.translate(molecule_center-average_coordinates)
    # 将分子添加到表面上
    system = surfacecopy + mol
    return system
'''random place'''
def random_place(size):
    x=size[0]
    y=size[1]
    #z=size[2]
    # 随机参数范围！！！
    x_range = [0,2.7*x]#Ru-Ru = 2.7A
    y_range = [0,2.7*y*((3**0.5)/2)]
    z_range = [1,5]#can change
    x_sv = np.random.uniform(x_range[0], x_range[1])
    y_sv = np.random.uniform(y_range[0], y_range[1])
    z_sv = np.random.uniform(z_range[0], z_range[1])
    sv = [x_sv,y_sv,z_sv]
    theta_z = np.degrees(np.random.uniform(0, 2*np.pi))
    varphi_y = np.degrees(np.random.uniform(0, 2*np.pi))
    return sv,theta_z,varphi_y
## 检查原子之间距离（H：0.5埃；other atom：1埃）
def check_dist_between_atoms(structure):
    cutoff_H = 0.9  # H-other atoms(including H)
    cutoff_other = 1.0 #other atoms - other atoms(both except H)
    for i in range(len(structure)):
        for j in range(i+1, len(structure)):
            atom_i = structure[i]
            atom_j = structure[j]
            element_i = atom_i.symbol
            element_j = atom_j.symbol
            dist = structure.get_distances(i,j,mic=True)
            if element_i == 'H' or element_j == 'H':
                if dist < cutoff_H:
                    print(f'原子 {element_i}:{i} 和原子 {element_j}:{j} 之间的距离为 {dist}埃, 小于截断值 {cutoff_H}')
                    return False
                else:
                    pass
            else:
                if dist < cutoff_other:
                    print(f'原子 {element_i}:{i} 和原子 {element_j}:{j} 之间的距离为 {dist}埃, 小于截断值 {cutoff_other}')
                    return False
                else:
                    pass
    return True

## 检查分子是否位于表面以上 & 吸附原子是否位于分子最下方
def check_molecule_over_surface(surface,mol,sv):
    z_max = max(surface.positions[:,2])+0.5#高于表面0.5A的距离
    #molecule_center = sv + np.array([0,0,z_max])
    z_min_mol = min(mol.positions[:,2])
    if z_min_mol < z_max:
        print(f'部分原子距离催化剂表面不到0.5埃')
        return False
    '''elif molecule_center[2] > z_min_mol:
        print(f'吸附原子未位于最靠近表面位置')
        return False'''
    return True
            
## 读取分子
def txt_to_dict(filename):
    dictionary = {}
    with open(filename, 'r') as file:
        for line in file:
            # 移除行首行尾的空白字符
            line = line.strip()
            # 忽略空行和注释行（假设注释行以'#'开头）
            if line and not line.startswith('#'):
                # 分割键和值，最多分割一次
                parts = line.split(':', 1)
                if len(parts) == 2:
                    key, value = parts
                    dictionary[key.strip()] = value.strip()
                else:
                    print(f"无法解析的行：{line}")
    return dictionary
## 建模
def build_random_system(element,size,moltxt,molfloder,save_path,random_mol_num):#save_path= mol_to_ad
    ru_surface = build_slab(element,size,10,2,save_path)
    adGroup_dict = txt_to_dict(moltxt)
    adGroup_namelist = list(adGroup_dict.keys())
    adGroup_filelist = list(adGroup_dict.values())
    group_num = len(adGroup_namelist)
    txt_name = save_path + '//floder_name.txt'
    with open(txt_name, 'w') as file:
        pass
    if 'ASEtoadG_output' in moltxt:
        for i in range(group_num):
            path_of_mol = molfloder+'//'+adGroup_filelist[i]
            addfile = list(adGroup_filelist[i].split('_')) 
            adGroup_mol = read(path_of_mol)
            adGroup_name = ast.literal_eval(adGroup_namelist[i])[-1]
            contains_star = '*' in adGroup_name
            if contains_star == True:
                adG_n = re.sub(r"\*", "-", adGroup_name)
            else:
                adG_n = re.sub(r"^", "--",adGroup_name)
            with open(txt_name, 'a') as file:
                species_file_floder_name = adG_n
                file.write(f'{species_file_floder_name}:{addfile[0]}\n')
            j = 1
            while j <= random_mol_num:#随机模型数量
                sv,theta_z,varphi_y = random_place(size)
                mol = rotate_mol(adGroup_mol,(theta_z,'z'),(varphi_y,'y'))
                system = place_mol_on_surface(mol,ru_surface,sv)
                if check_dist_between_atoms(system) == True and check_molecule_over_surface(ru_surface,mol,sv) == True:
                    floder_n = save_path+'/'+adG_n
                    file_n='POSCAR'
                    save_file(floder_n,file_n,system)
                    j=j+1
                else:
                    j=j
    elif 'SMItoASE_output' in moltxt:
        for i in range(group_num):
            path_of_mol = molfloder+'//'+adGroup_filelist[i]
            adGroup_mol = read(path_of_mol)
            adGroup_name = adGroup_namelist[i]
            with open(txt_name, 'a') as file:
                species_file_floder_name = adGroup_name
                file.write(f'{species_file_floder_name}:{adGroup_name}\n')
            j = 1
            while j <= random_mol_num:#随机模型数量
                sv,theta_z,varphi_y = random_place(size)
                mol = rotate_mol(adGroup_mol,(theta_z,'z'),(varphi_y,'y'))
                system = place_mol_on_surface(mol,ru_surface,sv)
                if check_dist_between_atoms(system) == True and check_molecule_over_surface(ru_surface,mol,sv) == True:
                    floder_n = save_path+'/'+adGroup_name+'/'+'struct_'+str(j)
                    #file_n=str(j)+'_'+adGroup_name+'.vasp'
                    file_n = 'POSCAR'
                    save_file(floder_n,file_n,system)
                    j=j+1
                else:
                    j=j



    