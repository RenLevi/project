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
function run_task() {
    local subdir=$1
    echo "Starting optimization in directory: $subdir"
    cd "$subdir"
    echo "Starting optimization" | tee -a resLog.out
    python mlp_calEnergy.py | tee -a resLog.out
    echo "" | tee -a resLog.out
    echo "Evaluation finished in directory: $dsubdir"
    echo "See resLog.out to check the results"
    cd ..
}
for dir in */;do
    echo "enter folder to optimization :$dir"
    cd "$dir"
    for subdir in struct_*/; do
        run_task "$subdir"
    done
    cd ..
    echo "finish optimization in folder :$dir"
done
echo "All jobs completed"

    
