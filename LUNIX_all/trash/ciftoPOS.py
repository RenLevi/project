import os
from ase.io import read, write

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
def cif2poscar(path1):
    # 定义文件夹路径
    source_folder = path1
    destination_folder = os.path.join(source_folder, 'POSCARs')

    # 创建 POSCARs 文件夹
    os.makedirs(destination_folder, exist_ok=True)

    # 初始化跳过的结构计数器
    skipped_count = 0

    # 遍历 source_folder 内的所有 CIF 文件
    for filename in os.listdir(source_folder):
        if filename.endswith('.cif'):
            cif_path = os.path.join(source_folder, filename)
            
            try:
                # 读取 CIF 文件
                structure = read(cif_path)
                
                # 生成 POSCAR 文件名
                poscar_filename = filename.replace('.cif', '.vasp')
                poscar_path = os.path.join(destination_folder, poscar_filename)
                
                # 将结构保存为 POSCAR 文件
                write(poscar_path, structure, format='vasp')
            except Exception as e:
                print(f"Skipped {filename}: {e}")
                skipped_count += 1

    print(f"Conversion complete. POSCAR files are saved in the POSCARs folder.")
    print(f"Total skipped structures: {skipped_count}")

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
        cif2poscar(target_directory)