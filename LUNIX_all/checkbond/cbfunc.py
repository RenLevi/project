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
    folder_dict = txt_to_dict(path1)
    namelist = list(folder_dict.keys())
    smileslist =list(folder_dict.values())
    txt1 = path2 +'checkbondpass.txt'
    txt2 = path2 +'wrong.txt'
    with open(txt1, 'w') as file:
        pass
    with open(txt2, 'w') as file:
        pass
    for foldername in namelist:
        fp = path2+foldername
        test = checkBonds()
        test.input(fp+'/opt.vasp')
        if test.CheckPBC() == True:
            test.AddAtoms()
            test.CheckAllBonds()
        else:
            pass
        output = BuildMol2Smiles(test)
        output.build()
        if output.smiles == foldername: #and output.ads != []:
            with open(txt1, 'a') as file:
                file.write(f'{foldername}:[{output.smiles},{output.smiles == foldername},{output.ads != []}]\n')
                print(f'{foldername} check pass')
        else:
            with open(txt2,'a') as file:
                file.write(f'{foldername}:[{output.smiles},{output.smiles == foldername},{output.ads != []}]\n')
                if output.smiles != foldername:
                    print(f'{foldername} have wrong with SMILES')
                if output.ads == []:
                    print(f'{foldername} have wrong with adsorption')


            
         