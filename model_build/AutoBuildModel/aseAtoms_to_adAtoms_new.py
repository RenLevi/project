'''v1.1.0'''
#input:txt of SMILES ===> output:.xyz of ads_molecules(Product,Reactant,Intermediate)
from openbabel import pybel
import numpy as np
from ase import Atoms
from ase.io import write
from rdkit import Chem
import os
import re
from SMILES_to_aseAtoms import *
#############################向输入的SMILES表达式中添加虚拟原子###################################################
def Bond_to_adGroup(smiles,adatomIdx):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    mw = Chem.RWMol(m)
    ghost_atom_index = mw.AddAtom(Chem.Atom(0))
    mw.AddBond(adatomIdx, ghost_atom_index, Chem.BondType.SINGLE)
    m_edit = mw.GetMol()
    Chem.SanitizeMol(m_edit)
    new_ad_smiles = Chem.MolToSmiles(m_edit)
    new_ad_smiles = add_brackets_around_letters(new_ad_smiles)#
    return new_ad_smiles 
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
        if bond > 1:
            if warp(atom.GetAtomicNum(),bond) == True:
                mayad.append(atom.GetIdx())
        elif bond < 1:
            print('Wrong with SMILES')
            break
        else:
            pass
    return mayad
def SMILES_to_mol_geo(smi):#SMILES格式正确
    # 使用Openbabel创建分子对象
    smiles = enumerate_smiles(smi)
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
    return smiles , molecule_geo
## 输出吸附基团的xyz文件
def SMILES_to_adGroup(smi,adAtom,file_format,save_path,txt_name):
    smiles = Bond_to_adGroup(smi[0],adAtom)
    smiles = enumerate_smiles(smiles)
    file_smi=re.sub(r"\*", "-", smiles)
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
    del molecule_geo[0]
    zeropointshift = molecule_geo.positions[0]
    molecule_geo.translate(-zeropointshift)
    file_name=smi[1]+'_'+file_smi+'.'+file_format#可更改文件保存路径以及格式
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    write(os.path.join(save_path, file_name), molecule_geo)
    with open(txt_name, 'a') as file:
        file.write(f'{[smi[1],smiles]}:{file_name}\n')
def SMILES_to_adMol(smi,file_format,save_path,txt_name):
        # 使用Openbabel创建分子对象
    smiles = enumerate_smiles(smi[0])
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
    file_name = smi[1]+'.'+file_format#可更改文件保存路径以及格式
    write(os.path.join(save_path, file_name), molecule_geo)
    with open(txt_name, 'a') as file:
        file.write(f'{[smi[1],smi[1]]}:{file_name}\n')
def adGroup_modeling(smiles_fromCN,file_format,save_path,txt_name):
    smiles, mol_geo = SMILES_to_mol_geo(smiles_fromCN[0])#得到规范化的smiles:[CH]C,和自由基模型
    smi =(smiles,smiles_fromCN[1])#([CH2],[H]C[H])
    adslist = find_ads_atoms(smi[0])
    if len(adslist) == 0:
        SMILES_to_adMol(smi,file_format,save_path,txt_name)
    else:
        for adatom in adslist:
            SMILES_to_adGroup(smi,adatom,file_format,save_path,txt_name)
    '''elif len(adslist) == 1:
        adatom = adslist[0]
        SMILES_to_adGroup(smi,adatom,file_format,save_path,txt_name)
    else:
        for i in range(len(adslist)):
            adatom = adslist[i]
            SMILES_to_adGroup(smi,adatom,file_format,save_path,txt_name)'''
    

def mollist_to_adGroup_files(mollist,file_format,save_path):#mollist=='mollist.txt
    
    """"
    
    DO NOT name your file as 'mollist.txt',if your file is NOT created by CNnet
    
    """

    if mollist == 'mollist.txt':
        mol_list = CatNetmol_to_SMILES(mollist)#展开的smiles (std_mol,CNmol)：([H][C][H],[H]C[H])
    else:
        mol_list = File_to_SMILES(mollist)
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