from aseAtoms_to_adAtoms_new import mollist_to_adGroup_files
from SMILES_to_aseAtoms import mollist_to_files
from molecule_ad_Cat import build_random_system
""""
    
    DO NOT name your file as 'mollist.txt',if your file is NOT created by CNnet
    
"""

'''生成吸附基团'''
txtName = 'mollist.txt'
fileFormat ='xyz'
savePath = "cal//output//"
element = 'Ru'
folderName = ['ASEtoadG_output','SMItoASE_output','mol_to_ad']

#mollist_to_adGroup_files(txtName,fileFormat,savePath+folderName[0])
#print(f"finish saving all files of adsorption Group in folder called {folderName[0]}")

'''生成自由基'''
mollist_to_files(txtName,fileFormat,savePath+folderName[1])
print(f"finish saving all files of free radical in folder called {folderName[1]}")

'''生成随机结构'''
build_random_system(element,(4,4,4),savePath+folderName[1]+'//species_name.txt',savePath+folderName[1]+'//species',savePath+folderName[2],50)
print(f"finish saving all files of adsorption Systems in folder called {folderName[2]}")