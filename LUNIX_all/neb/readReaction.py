from rdkit import Chem
from rdkit.Chem import rdmolops
from CheckNN import *
from ase.io import write
import numpy as np
import copy
from ase.data import covalent_radii, atomic_numbers
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
def angle_between_vectors(v1, v2):
    """使用NumPy的线性代数函数"""
    # 归一化向量
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    # 计算夹角的余弦值
    cos_theta = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
    return np.degrees(np.arccos(cos_theta))
def add_brackets_around_letters(cnmol:str):# 使用正则表达式替换不在[]中的字母字符，前后添加[]:example:[H]CO==>[H][C][O]
    result = re.sub(r'(?<!\[)([a-zA-Z])(?!\])', r'[\g<1>]', cnmol)
    return result
def subHH(STR):
    result = re.sub(r'\[HH\]', '[H]', STR)
    return result
def str2list(reaction:str):
    r1 = reaction.split(">")
    r2 = []
    for i in r1:
        i1 = i.split()
        r2.append(i1)
    return r2
def checkbond(reaction:list,bms1,bms2):
    mol1 = bms1.mol
    mol2 =bms2.mol
    reactiontype = reaction[1][0]
    addatom = reaction[1][-3]
    bondedatom = reaction[1][-1]
    def COMBINE(mol,add):
        # 创建分子时禁用化合价检查
        params = Chem.SmilesParserParams()
        params.removeHs = False  # 不自动移除氢原子
        params.sanitize = False  # 禁用所有检查
        addmol = Chem.MolFromSmiles(add_brackets_around_letters(add),params=params)
        combined_mol = rdmolops.CombineMols(mol, addmol)
        return combined_mol
    def addATOM():
        if reactiontype == 'Add':
            if len(addatom) > 1 and 'C'in addatom:
                return 'C'
            else:
                return addatom
        else:
            if addatom == 'O/OH':
                return 'O'
            else:
                return addatom
    def warp(cs12,add=addatom):
        check_mol = copy.deepcopy(cs12[1])
        check_mol = COMBINE(check_mol,add)
        bonds = cs12[0].GetBonds()
        AA=addATOM()
        aset = {AA,bondedatom}
        for bond in bonds:
            mol = copy.deepcopy(cs12[0])
            begin_atom_id = bond.GetBeginAtomIdx()
            end_atom_id = bond.GetEndAtomIdx()
            begin_atom = mol.GetAtomWithIdx(begin_atom_id)
            end_atom = mol.GetAtomWithIdx(end_atom_id)
            qset = {begin_atom.GetSymbol(),end_atom.GetSymbol()}
            if qset == aset:
                mol.RemoveBond(begin_atom_id, end_atom_id)
                if subHH(Chem.MolToSmiles(mol)) == Chem.MolToSmiles(check_mol):
                    return begin_atom.GetIdx(),end_atom.GetIdx(),Chem.MolToSmiles(cs12[0])
                else:
                    pass
            else:
                pass
        return False, False
    if reactiontype == 'Add':
        cs12 = (mol2,mol1)
        return warp(cs12)
    elif reactiontype == 'Remove':
        cs12 = (mol1,mol2)
        if addatom == 'O/OH':
            o1,o2 = warp(cs12,'O')
            if o1 == False and o2 == False:
                return warp(cs12,'OH')
            else:
                return o1,o2
        else:
            return warp(cs12)
            
def check_molecule_over_surface_and_not_cross_pbc(atoms):
    z_max = 17.415
    z_plist=[]
    nonmetals = ['H', 'He', 
                 'B', 'C', 'N', 'O', 'F', 'Ne',
                 'Si', 'P', 'S', 'Cl', 'Ar',
                 'Ge', 'As', 'Se', 'Br', 'Kr',
                 'Sb', 'Te', 'I', 'Xe',
                 'Po', 'At', 'Rn']
    for atom in atoms:
        symbol = atom.symbol
        if symbol in nonmetals:
            z_pos = atom.position[2]  # z坐标是position数组的第三个元素
            z_plist.append(z_pos)
    z_min = min(z_plist)
    if z_min < z_max:
        print(f'部分原子位于催化剂表面以下')
        return False
    else:    return True

def adjust_distance(readatoms, index1, index2,idlist,new_distance=0,delta=0,noads=False):
    """
    调整两个原子之间的距离
    
    参数:
        atoms: ASE Atoms 对象
        index1: 第一个原子的索引
        index2: 第二个原子的索引
        new_distance: 新的距离 (Å)
    """
    atoms = copy.deepcopy(readatoms)
    pos1 = atoms.positions[index1]
    pos2 = atoms.positions[index2]
    n1 = atoms.get_masses()[index1]
    n2 = atoms.get_masses()[index2]
    if noads == False:pass
    else:
        molIdxlist=[]
        for atom in atoms:
            if atom.symbol in ['C','H','O']:
                molIdxlist.append(atom.index)
            else:pass
        group = atoms[molIdxlist]
        v_important = pos2-pos1
        z= np.array([0,0,-1])
        theta = angle_between_vectors(v_important,z)
        axis_vz= np.cross(v_important,z)
        group.rotate(v=axis_vz,a=theta,center=pos1)
        val=v_important
        if are_vectors_parallel(val,np.array([0,0,-1])) == False:
            group.rotate(v=axis_vz,a=-2*theta,center=pos1)
        group.translate((0,0,17.5-pos1[2]))
        atoms.positions[molIdxlist] = group.positions
        print(f'{pos2-pos1}')
    r_1 = covalent_radii[int(n1)]
    r_2 = covalent_radii[int(n2)]
    new_distance = 1.5*(r_1 + r_2)
    vector = pos2 - pos1
    unit_vector = vector / np.linalg.norm(vector)
    v = unit_vector * new_distance
    pos2_new = v+pos1
    z1= 17+delta#pos1[2]+delta
    h = np.array([0,0,z1-pos2_new[2]])
    v13 = (v+h)*np.linalg.norm(v)/np.linalg.norm(v+h)
    v_final = copy.deepcopy(v13)
    # 移动第二个原子到新位置
    for id in idlist:
        atoms.positions[id] = atoms.positions[id]+v_final
    if noads == False:pass
    else:
        addgroup = atoms[idlist]
        v_important = pos2-pos1
        z= np.array([0,0,-1])
        theta = angle_between_vectors(v_important,z)
        axis_vz= np.cross(v_important,z)
        addgroup.rotate(v=axis_vz,a=theta,center=pos2)
        val=pos2-pos1
        if are_vectors_parallel(val,np.array([0,0,-1])) == False:
            addgroup.rotate(v=axis_vz,a=-2*theta,center=pos1)
        atoms.positions[idlist]=addgroup.positions
    return atoms
def check_neighbor(id,cb):
    idlist = []
    centeratom = cb.atoms[id]
    bonddict = centeratom.bonddict
    for atom in bonddict:
        idlist.append(atom.id)
    idlist.append(centeratom.id)
    return idlist
class readreaction():
    def __init__(self,file1,file2,reaction,noads=False):# file1> reaction > file2
        self.mol1 = file1
        self.mol2 = file2
        self.r = str2list(reaction)
        self.noads = noads
    def readfile(self):
        def warp(id,cb):
            ba = cb.atoms[id]
            if ba.elesymbol == bondedatom:
                return True
            else:return False
        CB1 = checkBonds()
        CB1.input(self.mol1)
        CB1.AddAtoms()
        CB1.CheckAllBonds()
        CB2 = checkBonds()
        CB2.input(self.mol2)
        CB2.AddAtoms()
        CB2.CheckAllBonds()
        BMS1 = BuildMol2Smiles(CB1)
        BMS1.build()
        BMS2 = BuildMol2Smiles(CB2)
        BMS2.build()
        begin_id,end_id,smilesFORcheck = checkbond(self.r,BMS1,BMS2)
        Bid_infile = begin_id +BMS1.metal 
        Eid_infile = end_id +BMS1.metal
        reactiontype = self.r[1][0]
        addatom = self.r[1][1]
        bondedatom = self.r[1][-1]
        if reactiontype == 'Add':
            CB = CB2
            self.nebFS = CB2.poscar
        else:
            CB = CB1
            self.nebFS = CB1.poscar
        if CB.adsorption == []:
            noads = True
        else:
            noads = False
        if warp(Bid_infile,CB) == True:
            idlist = check_neighbor(Eid_infile,CB)
            idlist.remove(Bid_infile)
            newmol = adjust_distance(CB.poscar,Bid_infile,Eid_infile,idlist,noads=noads)
            if check_molecule_over_surface_and_not_cross_pbc(newmol) == False:
                newmol = adjust_distance(CB.poscar,Bid_infile,Eid_infile,idlist,delta=1,noads=noads)
        else:
            idlist = check_neighbor(Bid_infile,CB)
            idlist.remove(Eid_infile)
            newmol = adjust_distance(CB.poscar,Eid_infile,Bid_infile,idlist,noads=noads)
            if check_molecule_over_surface_and_not_cross_pbc(newmol) == False:
                newmol = adjust_distance(CB.poscar,Eid_infile,Bid_infile,idlist,delta=1,noads=noads)
        self.nebIS = newmol
        self.check =smilesFORcheck 
    def save(self,path,format):
        # 保存为POSCAR文件（VASP格式）
        if format=='poscar' or 'POSCAR' or 'vasp':
            write(path+'IS.vasp', self.nebIS, format='vasp', vasp5=True)  # vasp5=True添加元素名称
            write(path+'FS.vasp', self.nebFS, format='vasp', vasp5=True)  # vasp5=True添加元素名称
        else:
            print('format should be .vasp')
        #test IS whether the bond ia breaked     
        test = checkBonds()
        test.input(path+'IS.vasp')
        if test.CheckPBC() == True:
            test.AddAtoms()
            test.CheckAllBonds()
        else:
            print('something wrong with pbc')
            #return TypeError   
        output = BuildMol2Smiles(test)
        output.build()
        if output.smiles == self.check:
            print(f'IS:{output.smiles},Check:{self.check}',output.smiles == self.check,'\n'
              'Error:the bond that should be breaked is not breaked'
              )
            #return ValueError
        if output.ads == []:
            print('Warning:there is not adsportion in model')        
            



        
if (__name__ == "__main__"):
        # 使用示例
    file1 = "cal/output/mol_to_ad/CO/opt.vasp"
    file2 = "cal/output/mol_to_ad/OCO/opt.vasp"
    reaction ='CO > Add O on C > OCO'
    PATH =  "test/"
    RR = readreaction(file1,file2,reaction)
    RR.readfile()
    RR.save(PATH,'POSCAR')






        

                






    

