import os
import glob 
from CheckNN import *
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
def find_files(directory, pattern):
    # 使用glob模块来查找匹配特定模式的文件
    files = glob.glob(os.path.join(directory, f'*{pattern}*'))
    return files


if __name__ == "__main__":
    # 示例用法
    path1 = 'cal/output/mol_to_ad/floder_name.txt'
    path2 = 'cal/output/mol_to_ad/'
    #path1 = 'cal/val/mol_to_ad/floder_name.txt'
    #path2 = 'cal/val/mol_to_ad/'
    folder_dict = txt_to_dict(path1)
    namelist = list(folder_dict.keys())
    smileslist =list(folder_dict.values())
    txt1 = path2 +'checkbondpass.txt'
    txt2 = path2 +'wrong.txt'
    with open(txt1, 'w') as file:
        pass
    with open(txt2, 'w') as file:
        pass
    for folder in range(len(namelist)):
        count_spilt = 0
        count_noads =0
        count_Hads = 0
        foldername = namelist[folder]
        sminame = smileslist[folder]
        for i in range(1,3):#
            fp = path2+foldername+'/'+'struct_'+str(i)
            test = checkBonds()
            test.input(fp+'/opt.vasp')
            if test.CheckPBC() == True:
                test.AddAtoms()
                test.CheckAllBonds()
            else:
                pass
            output = BuildMol2Smiles(test)
            output.build()
            set1=set()
            for a in output.ads:
                set1.add(a.elesymbol)
            if output.smiles == foldername  and output.ads != [] and set1 !={'H'}:
                with open(txt1, 'a') as file:
                    file.write(f'{foldername}:["{i}","{output.smiles}",{output.smiles == sminame},{set1},{set1!={'H'} and output.ads != []}]\n')
                    print(f'{foldername}_{i} check pass')
                    #break#
            else:
                if output.smiles != foldername:
                    count_spilt +=1
                if output.ads == []:
                    count_noads +=1
                if set1=={'H'}:
                    count_Hads +=1

        with open(txt2,'a') as file:
            f =2#
            if count_spilt == f:
                file.write(f'{foldername}:Dissociation\n')
                print(f'{foldername} have wrong with bonds')
            if count_noads == f:
                file.write(f'{foldername}:Not adsorbed\n')
                print(f'{foldername}_{i} have wrong with adsorption')
            if count_Hads == f:
                file.write(f'{foldername}:Not adsorbed|H ads\n')
                print(f'{foldername}_{i} have wrong with adsorption|H ads')
            
            if count_spilt!=f and count_noads!=f and count_Hads!=f:
                if count_noads+count_spilt==f:
                    file.write(f'{foldername}:noads & dissociation\n')
                    print(f'{foldername}_{i} have wrong with adsorption and dissociation')
                elif count_noads+count_Hads==f:
                    file.write(f'{foldername}:noads & H ads\n')
                    print(f'{foldername}_{i} have wrong with adsorption')
                elif count_Hads+count_spilt==f:
                    file.write(f'{foldername}:dissociation & H ads\n')
                    print(f'{foldername}_{i} have wrong with adsorption and dissociation')


            
         