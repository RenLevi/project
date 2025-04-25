#!/bin/bash
#SBATCH -J sequential_opt       # 作业名
#SBATCH -p wzhcnormal           # 队列名
#SBATCH -N 1                    # 节点数
#SBATCH --ntasks-per-node=28    # 每节点进程数
#SBATCH --cpus-per-task=1       # 每进程占用核心数
#SBATCH -o %j.out               # 标准输出文件
#SBATCH -e %j.err               # 错误输出文件

source ~/.bashrc
conda activate op

# 定义一个函数来执行任务
run_task() {
    local dir=$1
    echo "Starting optimization in directory: $dir"
    cd "$dir"
    echo "Starting optimization" | tee -a resLog.out
    python mlp_calEnergy.py | tee -a resLog.out
    echo "" | tee -a resLog.out
    echo "Evaluation finished in directory: $dir"
    echo "See resLog.out to check the results"
    cd ..
}

find "$parent_folder" -mindepth 1 -maxdepth 1 -type d | while read -r subfolder; do
    # 进入子文件夹
    cd "$subfolder" || exit
        # 遍历每个结构目录并顺序执行任务
    for dir in struct_*/; do
        run_task "$dir"
    done
    echo "All jobs completed"
    
    # 返回母文件夹
    cd "$parent_folder" || exit
done

