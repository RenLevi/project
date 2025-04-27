# -*- coding: utf-8 -*-
from SMILES_to_aseAtoms import mollist_to_files
from molecule_ad_Cat import build_random_system
import os
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
"""" 
    DO NOT name your file as 'mollist.txt',if your file is NOT created by CNnet
    
"""

'''生成吸附基团'''
main_folder = 'cal/output/'
txtName = 'wrong.txt'
fileFormat ='xyz'
savePath = "cal/val/"
element = 'Ru'
folderName = ['ASEtoadG_output/','SMItoASE_output/','mol_to_ad/']
folder_dict = txt_to_dict(main_folder+folderName[-1]+txtName)
smilist = list(folder_dict.keys())
if not os.path.exists(savePath):
        os.makedirs(savePath)
with open(savePath+'mollist.txt', 'w') as file:
    file.write(f'We total got {len(smilist)} mol that are not adsorbed or Dissociation:\n') 
for smi in smilist:
    with open(savePath+'mollist.txt', 'a') as file:
        file.write(f'{smi}\n')
#mollist_to_adGroup_files(txtName,fileFormat,savePath+folderName[0])
#print(f"finish saving all files of adsorption Group in folder called {folderName[0]}")
'''生成自由基'''
mollist_to_files(savePath+'mollist.txt',fileFormat,savePath+folderName[1])
print(f"finish saving all files of free radical in folder called {folderName[1]}")

'''生成随机结构'''
build_random_system(element,(4,4,4),savePath+folderName[1]+'species_name.txt',savePath+folderName[1]+'species',savePath+folderName[2],20)
print(f"finish saving all files of adsorption Systems in folder called {folderName[2]}")