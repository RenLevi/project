import subprocess
import os
import shutil
import glob

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
def copyFolder(source_folder,dest_folder):
    for item in os.listdir(source_folder):
        source_item = os.path.join(source_folder, item)
        dest_item = os.path.join(dest_folder, item)
        shutil.copy2(source_item, dest_item)

    print(f"文件夹已成功复制到 {dest_folder}")
def get_filename_from_path(filepath):
    """从文件路径中提取文件名（包含扩展名）"""
    return os.path.basename(filepath)
def find_files(directory, pattern):
    # 使用glob模块来查找匹配特定模式的文件
    files = glob.glob(os.path.join(directory, f'*{pattern}*'))
    return files
def read_file_line_by_line(file_path):#逐行读取txt文件并返回list数据
    mol_list=[]
    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            string = line.strip()  
            mol_list.append(string)
        mol_list.pop(0)
    return mol_list
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
def run_command_in_directory(directory, command):
    """
    在指定目录下运行命令
    
    参数:
        directory (str): 要切换到的目录路径
        command (str or list): 要执行的命令(字符串或列表形式)
    """
    # 保存当前工作目录
    original_dir = os.getcwd()
    
    try:
        # 切换到目标目录
        os.chdir(directory)
        print(f"当前工作目录已切换到: {os.getcwd()}")
        
        # 执行命令
        print(f"正在执行命令: {command}")
        result = subprocess.run(command, shell=isinstance(command, str), 
                               check=True, text=True, 
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # 输出命令执行结果
        print("命令输出:")
        print(result.stdout)
        if result.stderr:
            print("错误信息:")
            print(result.stderr)
            
    except subprocess.CalledProcessError as e:
        print(f"命令执行失败，返回码: {e.returncode}")
        print(f"错误信息: {e.stderr}")
    except Exception as e:
        print(f"发生错误: {str(e)}")
    finally:
        # 无论成功与否，都恢复原始工作目录
        os.chdir(original_dir)
        print(f"已恢复工作目录到: {original_dir}")


if __name__ == "__main__":
    path1 = '/public/home/ac877eihwp/renyq/cal/output/mol_to_ad/floder_name.txt'
    path2 = '/public/home/ac877eihwp/renyq/cal/output/mol_to_ad/'
    folder_dict = txt_to_dict(path1)
    namelist = list(folder_dict.keys())
    smileslist =list(folder_dict.values())
    for foldername in namelist:
        fp = path2+foldername
        copyFiles('/public/home/ac877eihwp/renyq/cal/jobsub.sh',fp)
        copyFiles('/public/home/ac877eihwp/renyq/cal/mlp_calEnergy.py',fp)

if __name__ == "__main__":
    # 示例用法
    path1 = '/public/home/ac877eihwp/renyq/cal/output/mol_to_ad/floder_name.txt'
    path2 = '/public/home/ac877eihwp/renyq/cal/output/mol_to_ad/'
    folder_dict = txt_to_dict(path1)
    namelist = list(folder_dict.keys())
    smileslist =list(folder_dict.values())
    for foldername in namelist:
        fp = path2+foldername
        target_directory = fp  # 替换为你的目标目录
        command_to_run = "qsub jobsub.sh"  # 替换为你想运行的命令
    
    # 也可以使用列表形式传递命令(推荐，更安全)
    # command_to_run = ["ls", "-l"]
        run_command_in_directory(target_directory, command_to_run)
        
















