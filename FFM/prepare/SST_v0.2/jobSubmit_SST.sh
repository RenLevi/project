#!/bin/bash
#SBATCH -J vasp544  #作业名
#SBATCH -p ihicnormal #队列名
#SBATCH -N 1 #节点数
#SBATCH --ntasks-per-node=1    #节点进程数
#SBATCH --cpus-per-task=28
#SBATCH --exclusive  ##独占节点，按节点计费
#SBATCH -o %j.out
#SBATCH -e %j.err


cd /public/home/ac877eihwp/renyq/TEST/cleanData/4_batch
source /public/home/ac877eihwp/.bashrc
conda activate nequip

python main.py merged.xyz
