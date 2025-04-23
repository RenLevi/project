import MLPopt_suppotfunc as MSF
import os
pathtxt = '/work/home/ac877eihwp/renyq/cal/output/mol_to_ad/floder_name.txt'
pathads = '/work/home/ac877eihwp/renyq/cal/output/mol_to_ad/'
path_of_support = ['/work/home/ac877eihwp/renyq/LUNIX_all/mlp_opt/mlp_calEnergy.py','/work/home/ac877eihwp/renyq/LUNIX_all/mlp_opt/jobsubopt.sh']
folder_dict = MSF.txt_to_dict(pathtxt)
namelist = list(folder_dict.keys())
smileslist =list(folder_dict.values())
for foldername in namelist:
    for i in range(1,21):
        fp =pathads+foldername+'/'+str(i)
        MSF.copyFiles(path_of_support[0],fp)
        MSF.copyFiles(path_of_support[1],fp)

for foldername in namelist:
    for i in range(1,21):
        fp = pathads+foldername+'/'+str(i)
        target_directory = fp  # 替换为你的目标目录
        command_to_run = "sbatch jobsubopt.sh"  # 替换为你想运行的命令
    
    # 也可以使用列表形式传递命令(推荐，更安全)
    # command_to_run = ["ls", "-l"]
        MSF.run_command_in_directory(target_directory, command_to_run)


    

