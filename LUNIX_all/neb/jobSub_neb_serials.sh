#!/bin/bash
#SBATCH -J sequential_neb       # 作业名
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
    local dir=$1
    echo "Starting optimization in directory: $dir"
    echo "Starting optimization" | tee -a resLog.out
    python test_neb.py | tee -a resLog.out
    echo "" | tee -a resLog.out
    echo "Evaluation finished"
    echo "See resLog.out to check the results"
}
for dir in */;do
    echo "enter folder to start neb to search TS :$dir"
    cd "$dir"
    run_task "$dir"
    cd ..
    echo "finish searching TS in folder :$dir"
done
echo "All jobs completed"

    
