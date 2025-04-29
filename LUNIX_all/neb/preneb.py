import readReaction as rR
import os
import shutil
import ast
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
def read_file_line_by_line(file_path):#逐行读取txt文件并返回list数据
    reaction_list=[]
    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            string = line.strip()  
            reaction_list.append(string)
        reaction_list.pop(0)
    return reaction_list
def copyFiles(source_file,dest_folder):
# 源文件路径source_file = '/path/to/source/file.txt'
# 目标文件夹路径dest_folder = '/path/to/destination/folder'
# 确保目标文件夹存在，如果不存在则创建
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)
    # 目标文件路径
    dest_file = os.path.join(dest_folder, os.path.basename(source_file))
    # 复制文件
    try:
        shutil.copy2(source_file, dest_file)
        print(f"文件已成功复制到 {dest_file}")
    except IOError as e:
        print(f"无法复制文件. {e}")
    except Exception as e:
        print(f"发生错误: {e}")

path0 = 'reactionslist.txt'
path1 = 'cal/output/mol_to_ad/floder_name.txt'
path2_1 = 'cal/output/mol_to_ad/checkbondpass.txt'
path2_2 = 'cal/val/mol_to_ad/checkbondpass.txt'
path3 = 'cal/val/mol_to_ad/wrong.txt'
path4 = 'cal/'
path5 = ['cal/output/mol_to_ad/','cal/val/mol_to_ad/']
mainfolder = path4+'neb/'
os.makedirs(mainfolder, exist_ok=True)  # exist_ok=True 避免目录已存在时报错
folder_dict = txt_to_dict(path1)
dict_pass1 = txt_to_dict(path2_1)
dict_pass2 = txt_to_dict(path2_2)
wrongdict = txt_to_dict(path3)
wrong = list(wrongdict.keys())
reaction_list = read_file_line_by_line(path0)
with open(mainfolder+'foldername.txt', 'w') as file:
    pass
for reaction in reaction_list:
    rlist = rR.str2list(reaction)
    if rlist[0][0] in wrong or rlist[-1][0] in wrong:
        if rlist[0][0] in wrong:
            initial=rlist[0][0]
            state_ini=wrongdict[initial]
        else:
            state_ini='Pass'
        if rlist[-1][0] in wrong:
            final= rlist[-1][0]
            state_fin=wrongdict[final]
        else:
            state_fin='Pass'
        if  state_ini in ['Dissociation','wrong'] or state_fin in ['Dissociation','wrong']:
            pass
        else:#not ads
            subfolder = mainfolder + str(rlist[0][0])+'_'+str(rlist[-1][0])+'/'
            with open(mainfolder+'foldername.txt', 'a') as file:
                reaction_floder_name = str(rlist[0][0])+'_'+str(rlist[-1][0])
                file.write(f'{reaction_floder_name}:{reaction}\n')
            os.makedirs(subfolder, exist_ok=True)
            file1 = path5[0]+str(rlist[0][0])+'/struct_1'+'/opt.vasp'
            file2 = path5[0]+str(rlist[-1][0])+'/struct_1'+'/opt.vasp'
            RR = rR.readreaction(file1,file2,reaction,noads=True)
            RR.readfile()
            RR.save(subfolder,'POSCAR')
            copyFiles('LUNIX_all/neb/test_neb.py',subfolder)
            copyFiles('LUNIX_all/neb/jobsubneb.sh',subfolder)
    else:
        subfolder = mainfolder + str(rlist[0][0])+'_'+str(rlist[-1][0])+'/'
        with open(mainfolder+'foldername.txt', 'a') as file:
            reaction_floder_name = str(rlist[0][0])+'_'+str(rlist[-1][0])
            file.write(f'{reaction_floder_name}:{reaction}\n')
        os.makedirs(subfolder, exist_ok=True)
        if rlist[0][0] in dict_pass1:
            val1=0
            list1 =ast.literal_eval(dict_pass1[rlist[0][0]])
        else:
            val1=1
            list1 =ast.literal_eval(dict_pass2[rlist[0][0]])
        if rlist[-1][0] in dict_pass1:
            val2=0
            list2 =ast.literal_eval(dict_pass1[rlist[-1][0]])
        else:
            val2=1
            list2 =ast.literal_eval(dict_pass2[rlist[-1][0]])
        file1 = path5[val1]+str(rlist[0][0])+'/struct_'+list1[0]+'/opt.vasp'
        file2 = path5[val2]+str(rlist[-1][0])+'/struct_'+list2[0]+'/opt.vasp'
        RR = rR.readreaction(file1,file2,reaction)
        RR.readfile()
        RR.save(subfolder,'POSCAR')
        copyFiles('LUNIX_all/neb/test_neb.py',subfolder)
        copyFiles('LUNIX_all/neb/jobsubneb.sh',subfolder)
#copyFiles('LUNIX_all/neb/jobsubneb.sh',)





