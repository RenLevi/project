import subprocess
import os
import shutil
import glob

def move_specific_file(source_folder, destination_folder, file_name):
    """
    将源文件夹中的特定文件移动到目标文件夹。

    :param source_folder: 源文件夹路径
    :param destination_folder: 目标文件夹路径
    :param file_name: 要移动的文件名
    """
    # 检查源文件夹是否存在
    if not os.path.exists(source_folder):
        print(f"源文件夹 {source_folder} 不存在。")
        return False

    # 构建源文件和目标文件的完整路径
    source_file_path = os.path.join(source_folder, file_name)
    destination_file_path = os.path.join(destination_folder, file_name)

    # 检查文件是否存在
    if not os.path.isfile(source_file_path):
        print(f"文件 {source_file_path} 不存在。")
        return False

    # 使用shutil.move移动文件
    shutil.move(source_file_path, destination_file_path)
    print(f"文件 {source_file_path} 已移动到 {destination_file_path}")
    return True

def run_file_in_directory( target_dir,file_name,command):
    """
    进入指定目录，运行指定文件，然后返回到原始目录。

    :param current_dir: 当前工作目录
    :param target_dir: 目标目录
    :param file_name: 要运行的文件名
    """
    try:
        # 保存当前工作目录
        original_dir = os.getcwd()

        # 进入目标目录
        os.chdir(target_dir)
        print(f"已进入目录：{target_dir}")

        # 构建文件的完整路径
        full_file_path = os.path.join(target_dir, file_name)
        
        # 检查文件是否存在
        if not os.path.isfile(full_file_path):
            print(f"文件 {full_file_path} 不存在。")
            return

        # 运行文件
        print(f"正在运行命令：{command}")
        subprocess.run([command], check=True)

    finally:
        # 返回到原始目录
        os.chdir(original_dir)
        print(f"已返回到原始目录：{original_dir}")

def add_jobSubmit_serial_to_nth_batch(batch_num,source_folder):
    for i in range(1,batch_num+1):
        destination_folder = str(i)+'_batch'
        move_specific_file(source_folder,destination_folder,'jobSubmit_serial.sh')
def run_jobSubmit_serial_to_nth_batch(batch_num):
    for i in range(1,batch_num+1):
        destination_folder = str(i)+'_batch'
        run_file_in_directory(destination_folder,'jobSubmit_serial.sh','qsub jobSubmit_serial.sh')

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
        
















