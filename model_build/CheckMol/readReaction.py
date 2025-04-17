from rdkit import Chem
from rdkit.Chem import rdmolops
from CheckNN import *
from ase.io import write
import numpy as np
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
    smi1 = reaction[0][0]
    smi2 = reaction[-1][0]
    mol1 = bms1.mol
    mol2 =bms2.mol
    if {smi1,smi2}!={bms1.smiles,bms2.smiles}:
        return False,False
    reactiontype = reaction[1][0]
    addatom = reaction[1][1]
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
                    return begin_atom.GetIdx(),end_atom.GetIdx(),Chem.MolToSmiles(check_mol)
                else:
                    pass
            else:
                pass
        return False, False

    if reactiontype == 'Add':
        cs12 = (mol2,mol1)
        return warp(cs12)
    elif reactiontype == 'Reamove':
        cs12 = (mol1,mol2)
        if addatom == 'O/OH':
            o1,o2 = warp(cs12,'O')
            if o1 == False and o2 == False:
                return warp(cs12,'OH')
            else:
                return o1,o2
            
def check_molecule_over_surface_and_not_cross_pbc(atoms):
    z_max = 16.415
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

def adjust_distance(atoms, index1, index2,idlist,new_distance,delta=0):
    """
    调整两个原子之间的距离
    
    参数:
        atoms: ASE Atoms 对象
        index1: 第一个原子的索引
        index2: 第二个原子的索引
        new_distance: 新的距离 (Å)
    """
    pos1 = atoms.positions[index1]
    pos2 = atoms.positions[index2]
    
    vector = pos2 - pos1
    unit_vector = vector / np.linalg.norm(vector)
    v = unit_vector * new_distance
    pos2_new = v+pos1
    z1= pos1[2]+delta
    h = np.array([0,0,z1-pos2_new[2]])
    v13 = (v+h)*np.linalg.norm(v)/np.linalg.norm(v+h)
    v_final = copy.deepcopy(v13)
    newatoms = copy.deepcopy(atoms)
    # 移动第二个原子到新位置
    for id in idlist:
        newatoms.positions[id] = newatoms.positions[id]+v_final
    return newatoms
def check_neighbor(id,cb):
    idlist = []
    centeratom = cb.atoms[id]
    bonddict = centeratom.bonddict
    for atom in bonddict:
        idlist.append(atom.id)
    idlist.append(centeratom.id)
    return idlist
class readreaction():
    def __init__(self,file1,file2,reaction):# file1> reaction > file2
        self.mol1 = file1
        self.mol2 = file2
        self.r = str2list(reaction)
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
        print(begin_id,end_id)
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
        if warp(Bid_infile,CB) == True:
            idlist = check_neighbor(Eid_infile,CB)
            idlist.remove(Bid_infile)
            newmol = adjust_distance(CB.poscar,Bid_infile,Eid_infile,idlist,2)#2埃
            if check_molecule_over_surface_and_not_cross_pbc(newmol) == False:
                newmol = adjust_distance(CB.poscar,Bid_infile,Eid_infile,idlist,2,1)
        else:
            idlist = check_neighbor(Bid_infile,CB)
            idlist.remove(Eid_infile)
            newmol = adjust_distance(CB.poscar,Eid_infile,Bid_infile,idlist,2)
            if check_molecule_over_surface_and_not_cross_pbc(newmol) == False:
                newmol = adjust_distance(CB.poscar,Eid_infile,Bid_infile,idlist,2,1)
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
        output = BuildMol2Smiles(test)
        output.build()
        print(f'IS:{output.smiles}',output.smiles == self.check,'\n'
              'Error:the bond that should be breaked is not breaked'
              )
        if output.ads == []:
            print('Warning:there is not adsportion in model')        
            



        
if (__name__ == "__main__"):
        # 使用示例
    file1 = "C:/Users/renyq/Desktop/result/neb/Ru_[H]OC([H])O/nequipOpt.traj"
    file2 = 'C:/Users/renyq/Desktop/result/neb/Ru_[H]OC([H])([H])O/nequipOpt.traj'
    reaction ='[H]OC([H])O > Add H on C > [H]OC([H])([H])O'
    PATH =  "C:/Users/renyq/Desktop/result/nebs/"
    RR = readreaction(file1,file2,reaction)
    RR.readfile()
    RR.save(PATH,'POSCAR')






        

                






    

