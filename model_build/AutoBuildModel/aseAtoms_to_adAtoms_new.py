'''v1.1.0'''
#input:txt of SMILES ===> output:.xyz of ads_molecules(Product,Reactant,Intermediate)
from openbabel import pybel
import openbabel as ob
import numpy as np
from ase import Atoms
from ase.io import write
from rdkit import Chem
import os
import re
import copy
def are_vectors_parallel(v1, v2, tol=1e-6):
    """
    检查两向量是否方向相同（或相反）。
    返回:
        True  (方向相同: 点积 ≈ 1)
        True  (方向相反: 点积 ≈ -1)
        False (其他情况)
    """
    v1_unit = v1 / np.linalg.norm(v1)
    v2_unit = v2 / np.linalg.norm(v2)
    dot_product = np.dot(v1_unit, v2_unit)
    return np.isclose(abs(dot_product), 1.0, atol=tol)
def rotate_mol(mol,rotate_matrix1):#rotate molecule
    molcopy = copy.deepcopy(mol)
    (angle1,axis1) = rotate_matrix1
    molcopy.rotate(angle1,axis1,center=(0,0,0))
    return molcopy
def angle_between_vectors(v1, v2):
    """使用NumPy的线性代数函数"""
    # 归一化向量
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    # 计算夹角的余弦值
    cos_theta = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
    return np.degrees(np.arccos(cos_theta))
def enumerate_smiles(smiles):#列举同一分子smiles的等价表达,得到规范化的smiles
    mol = Chem.MolFromSmiles(smiles)
    smi = Chem.MolToSmiles(mol,doRandom=False,canonical=False)
    return smi
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
def File_to_SMILES(filename):
    molecule_list = read_file_line_by_line(filename)
    return molecule_list
def Bond_to_adGroup(smiles,adatomIdx):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    mw = Chem.RWMol(m)
    ghost_atom_index = mw.AddAtom(Chem.Atom(0))
    mw.AddBond(adatomIdx, ghost_atom_index, Chem.BondType.SINGLE)
    m_edit = mw.GetMol()
    Chem.SanitizeMol(m_edit)
    new_ad_smiles = Chem.MolToSmiles(m_edit)
    new_ad_smiles = enumerate_smiles(add_brackets_around_letters(new_ad_smiles))
    return new_ad_smiles 
def Bond_to_adMol_regadless_volence(smiles,adatomIdx):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    mw = Chem.RWMol(m)
    ghost_atom_index = mw.AddAtom(Chem.Atom(16))
    mw.AddBond(adatomIdx, ghost_atom_index, Chem.BondType.SINGLE)
    m_edit = mw.GetMol()
    m_edit.UpdatePropertyCache(strict=False)  # 不严格检查化合价
    Chem.SanitizeMol(m_edit, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)  # 跳过化合价检查
    # 输出 SMILES（可能不合理，但 RDKit 不会报错）
    smi = add_brackets_around_letters(Chem.MolToSmiles(m_edit,rootedAtAtom=ghost_atom_index))
    return smi #*[C]([H])([H])([H])[H]
## 查询吸附原子
def find_ads_atoms(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    mayad = []
    def warp(atomic_num,bondnum):
        periodic_table = Chem.GetPeriodicTable()
        possible_valence = periodic_table.GetDefaultValence(atomic_num)
        if bondnum < possible_valence:
            return True
        else:
            return False
    for atom in mol.GetAtoms():
        bond = atom.GetDegree()
        if warp(atom.GetAtomicNum(),bond) == True:
            mayad.append(atom.GetIdx())
        else:
            pass
    return mayad
def find_mol_ad_atoms(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    mayadNO = []
    mayadC = []
    def warp(atom):
        if atom.GetSymbol() == 'O' or atom.GetSymbol() == 'N':
            return 1
        elif atom.GetSymbol() == 'C':
            return 2
        elif atom.GetSymbol() == 'H':
            return 0
        else:
            return 3
    for atom in mol.GetAtoms():
        backvalue = warp(atom)
        if backvalue == 1:
            mayadNO.append(atom.GetIdx())
        elif backvalue == 2:
            mayadC.append(atom.GetIdx())
        else:
            pass
    if mayadNO != []:
        return mayadNO,'N&O'
    else:
        return [mayadC[0]],'C'
## 输出吸附基团的xyz文件
def SMILES_to_adGroup(smi,sy,adAtom,file_format,save_path,txt_name):
    if sy == None:
        smiles = Bond_to_adGroup(smi[0],adAtom)
    elif sy == 'C':
        smiles = Bond_to_adMol_regadless_volence(smi[0],adAtom)
    elif sy == 'N&O':
        smiles = Bond_to_adMol_regadless_volence(smi[0],adAtom)  
    if '*' in smiles:
        file_smi=re.sub(r"\*", "-", smiles)
    else:
        file_smi=re.sub(r"\[S\]", "--", smiles)
    molecule = pybel.readstring("smi", smiles)
    molecule.make3D(forcefield='mmff94', steps=100)
    # 创建ASE的Atoms对象
    molecule_geo = Atoms()
    #del molecule_geo[[molecule_geo.index for atom in molecule_geo if molecule_geo.symbol == 'X']]
    molecule_geo.set_cell(np.array([[6.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 6.0]]))
    molecule_geo.set_pbc((True, True, True))
    # 将Openbabel分子对象添加到ASE的Atoms对象中
    for atom in molecule:
        atom_type = atom.atomicnum
        atom_position = np.array([float(i) for i in atom.coords])
        molecule_geo.append(atom_type)
        molecule_geo.positions[-1] = atom_position

    v_important = molecule_geo.positions[0]-molecule_geo.positions[1]
    z= np.array([0,0,-1])
    theta = angle_between_vectors(v_important,z)
    axis_vz= np.cross(v_important,z)
    zeropointshift = molecule_geo.positions[1]
    molecule_geo.translate(-zeropointshift)
    molecule_geo = rotate_mol(molecule_geo,(theta,axis_vz))
    val=molecule_geo.positions[0]-molecule_geo.positions[1]
    if are_vectors_parallel(val,np.array([0,0,-1])) == False:
        molecule_geo = rotate_mol(molecule_geo,(-2*theta,axis_vz))
    del molecule_geo[0]
    file_name=smi[1]+'_'+file_smi+'.'+file_format#可更改文件保存路径以及格式
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    write(os.path.join(save_path, file_name), molecule_geo)
    with open(txt_name, 'a') as file:
        file.write(f'{[smi[1],re.sub(r"\[S\]", "", smiles)]}:{file_name}\n')
def adGroup_modeling(smiles_fromCN,file_format,save_path,txt_name):
    smi =(enumerate_smiles(add_brackets_around_letters(smiles_fromCN)),smiles_fromCN)#([CH2]O,[H]C[H]O[H])
    adslist = find_ads_atoms(smi[0])
    moladlist,sy = find_mol_ad_atoms(smi[0])
    if adslist == []:
        inputlist = moladlist
    else:
        inputlist = adslist
        sy = None
    for adatom in inputlist:
        SMILES_to_adGroup(smi,sy,adatom,file_format,save_path,txt_name)

def mollist_to_adGroup_files(mollist,file_format,save_path):#mollist=='mollist.txt
    mol_list = File_to_SMILES(mollist)#[H]OC([H])O
    species_saveFiles = save_path+'//species'
    if not os.path.exists(species_saveFiles):
        os.makedirs(species_saveFiles)
    txt_name = save_path + '//species_name.txt'
    with open(txt_name, 'w') as file:
        pass
    for smiles in mol_list:
        adGroup_modeling(smiles,file_format,species_saveFiles,txt_name)
'''Example:'''
#mollist_to_adGroup_files('mollist.txt','xyz','ASE//1.0.3//ASEtoadG_output')