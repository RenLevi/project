'''v1.1.0'''
#input:txt of SMILES ===> output:.xyz of molecules(Product,Reactant,Intermediate)
from openbabel import pybel
import numpy as np
from ase import Atoms
from ase.io import write
from rdkit import Chem
import re
import os

def check_input_SMILES(smiles):#检查输入的smiles是否符合标准
    # 尝试将SMILES字符串转换为分子对象
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False  # SMILES格式不正确
    return True  # SMILES格式正确

def enumerate_smiles(smiles):#列举同一分子smiles的等价表达,得到规范化的smiles
    mol = Chem.MolFromSmiles(smiles)
    smi = Chem.MolToSmiles(mol,doRandom=False,canonical=False)
    return smi

def SMILES_to_aseAtoms(smi_fromCN,file_format,save_path,txt_name):#SMILES格式需正确
    # 使用Openbabel创建分子对象
    smiles = enumerate_smiles(smi_fromCN[0])
    molecule = pybel.readstring("smi", smiles)
    molecule.make3D(forcefield='mmff94', steps=100)
    # 创建ASE的Atoms对象
    molecule_geo = Atoms()
    molecule_geo.set_cell(np.array([[6.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 6.0]]))
    molecule_geo.set_pbc((True, True, True))
    # 将Openbabel分子对象添加到ASE的Atoms对象中
    for atom in molecule:
        atom_type = atom.atomicnum
        atom_position = np.array([float(i) for i in atom.coords])
        molecule_geo.append(atom_type)
        molecule_geo.positions[-1] = atom_position
    # 调整分子位置,确定自由基位置
    free_radical_position = find_free_radical(smiles,molecule_geo)
    molecule_geo.translate(-free_radical_position)

    #save xyz file
    file_name = smi_fromCN[1]+'.'+file_format#可更改文件保存路径以及格式
    write(os.path.join(save_path, file_name), molecule_geo)
    with open(txt_name, 'a') as file:
        file.write(f'{smi_fromCN[1]}:{file_name}\n')

# 建模函数！！！最终步
def molecule_modeling(smiles,file_format,save_path,txt_name):#SMILES,保存格式，保存相对路径（相对于工作目录）
    if check_input_SMILES(smiles[0]) == False:
        print(f"The SMILES'{smiles[0]}'are incorrect")
    else:
        SMILES_to_aseAtoms(smiles,file_format,save_path,txt_name)
def mollist_to_files(mollist,file_format,save_path):#mollist=='mollist.txt
    
    """"
    
    DO NOT name your file as 'mollist.txt',if your file is NOT created by CNnet
    
    """

    if mollist == 'mollist.txt':
        mol_list = CatNetmol_to_SMILES(mollist)#[H][C][H]
    else:
        mol_list = File_to_SMILES(mollist)
    species_saveFiles = save_path+'//species'
    if not os.path.exists(species_saveFiles):
        os.makedirs(species_saveFiles)
    txt_name = save_path + '//species_name.txt'
    with open(txt_name, 'w') as file:
        pass
    for smiles in mol_list:
        molecule_modeling(smiles,file_format,species_saveFiles,txt_name)
# 判断自由基位点，非建模直接相关的辅助函数
def smiles_to_chemical_graph(smiles):
    # 从SMILES字符串创建分子对象
    mol = Chem.MolFromSmiles(smiles)
    # 添加氢原子
    mol = Chem.AddHs(mol)
    # 获取原子数量
    num_atoms = mol.GetNumAtoms()
    # 初始化原子种类列表
    atom_types = []
    # 初始化边列表，包含原子编号
    edges = []
    # 初始化邻接矩阵
    adjacency_matrix = np.zeros((num_atoms, num_atoms), dtype=int)
    # 遍历分子中的所有原子
    for i, atom in enumerate(mol.GetAtoms()):
        # 获取原子的符号（如'C'表示碳，'O'表示氧，'H'表示氢）
        atom_types.append((i, atom.GetSymbol()))
    # 遍历分子中的所有键
    for bond in mol.GetBonds():
        # 获取键的起始和结束原子索引
        begin_atom_idx = bond.GetBeginAtomIdx()
        end_atom_idx = bond.GetEndAtomIdx()
         # 将边添加到边列表中
        edges.append((begin_atom_idx, end_atom_idx))
        # 获取键的类型（单键、双键、三键或芳香键），并转换为整数
        bond_type = bond.GetBondType()
        # 在邻接矩阵中标记成键关系
        # 单键用1表示，双键用2表示，三键用3表示，芳香键用4表示
        bond_value = {
            Chem.BondType.SINGLE: 1,
            Chem.BondType.DOUBLE: 2,
            Chem.BondType.TRIPLE: 3,
            Chem.BondType.AROMATIC: 4
        }.get(bond_type, 1)  # 默认为单键
        adjacency_matrix[begin_atom_idx, end_atom_idx] = bond_value
        adjacency_matrix[end_atom_idx, begin_atom_idx] = bond_value  # 无向图，需要对称标记
    return atom_types, adjacency_matrix,edges#list,array,list

def count_nonzero_in_row(matrix):
    non_zero_counts = np.count_nonzero(matrix,axis=1)
    return non_zero_counts

def check_atoms_bonds(atom_types,bondarray):
    free_radical_position_list=[]
    for i in range(len(atom_types)):
        atom = atom_types[i]
        bonds = bondarray[i]
        if atom[1] == 'C' and bonds < 4:
            free_radical_position_list.append(atom)
        elif atom[1] == 'N' and bonds < 3:
            free_radical_position_list.append(atom)
        elif atom[1] == 'O' and bonds < 2:
            free_radical_position_list.append(atom)
        #elif atom[1] == 'H' and bonds < 1:
        #    free_radical_position_list.append(atom)
        else:
            pass
    return free_radical_position_list

def find_free_radical(smiles,molecule_geo):
    atom_types, adjacency_matrix,edges = smiles_to_chemical_graph(smiles)
    bondlist = count_nonzero_in_row(adjacency_matrix)
    freelist=check_atoms_bonds(atom_types,bondlist)
    if len(freelist) == 1:
        free_radical_position = molecule_geo.positions[freelist[0][0]]
    elif len(freelist) > 1:
        num = len(freelist)
        sum_position = molecule_geo.positions[freelist[0][0]]
        for i in range(1,len(freelist)):
            atom_ids = freelist[i][0]
            atom_position = molecule_geo.positions[atom_ids]
            sum_position = sum_position + atom_position
        free_radical_position = sum_position/num
    else:
        free_radical_position = molecule_geo.positions[0]
    return free_radical_position

# 修正非标准的SMILES表达式，非直接与建模相关的辅助函数
def add_brackets_around_letters(cnmol:str):# 使用正则表达式替换不在[]中的字母字符，前后添加[]:example:[H]CO==>[H][C][O]
    result = re.sub(r'(?<!\[)([a-zA-Z])(?!\])', r'[\g<1>]', cnmol)
    return result

def read_file_line_by_line(file_path):#逐行读取txt文件并返回list数据
    mol_list=[]
    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            string = line.strip()  
            mol_list.append(string)
        mol_list.pop(0)
    return mol_list
def CatNetmol_to_SMILES(filename):#将CatNet的输出文件中的非标准SMILES转为标准的SMILES 
    #默认原子之间成键均为单健,仍待考虑
    #得到分子数据的列表
    molecule_list = read_file_line_by_line(filename)
    std_mol_list = []
    for CNmol in molecule_list:
        std_mol = add_brackets_around_letters(CNmol)
        std_mol_list.append((std_mol,CNmol))#([H][C][H],[H]C[H])
    return std_mol_list
def File_to_SMILES(filename):
    molecule_list = read_file_line_by_line(filename)
    return molecule_list
'''
缺少读取mollist并输出文件的总函数
缺乏优化分子结构
旋转mode跟到后续放置分子于表面的步骤
make3D的forcefield
'xyz','./output/'
2024.11.5note:
\start
输入正确输出报错的问题无法解决
吸附问题的解决:dummy atom or ghost atom
\end
'''
########################################################
'''Example:'''
#mollist_to_files('mollist.txt','xyz','SMItoASE_output')
