import MLPopt_suppotfunc as MSF
import os
pathtxt = '/work/home/ac877eihwp/renyq/cal/val/mol_to_ad/floder_name.txt'
pathads = '/work/home/ac877eihwp/renyq/cal/val/mol_to_ad/'
#pathtxt = '/work/home/ac877eihwp/renyq/cal/output/mol_to_ad/floder_name.txt'
#pathads = '/work/home/ac877eihwp/renyq/cal/output/mol_to_ad/'
path_of_support = ['/work/home/ac877eihwp/renyq/LUNIX_all/mlp_opt/mlp_calEnergy.py','/work/home/ac877eihwp/renyq/LUNIX_all/mlp_opt/jobSubmit_serial.sh']
folder_dict = MSF.txt_to_dict(pathtxt)
namelist = list(folder_dict.keys())
smileslist =list(folder_dict.values())
for foldername in namelist:
    fp =pathads+foldername+'/'
    for i in range(1,21):#需改
        MSF.copyFiles(path_of_support[0],fp+'struct_'+str(i))
        
MSF.copyFiles(path_of_support[1],pathads)

#for foldername in namelist:
target_directory = pathads  # 替换为你的目标目录
command_to_run = "sbatch jobSubmit_serial.sh"  # 替换为你想运行的命令
    
    # 也可以使用列表形式传递命令(推荐，更安全)
    # command_to_run = ["ls", "-l"]
MSF.run_command_in_directory(target_directory, command_to_run)


    

