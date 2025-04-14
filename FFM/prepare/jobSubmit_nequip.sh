#!/bin/bash
#PBS -N nequip_revico
#PBS -l nodes=1:ppn=20
#PBS -l gpus=1
#PBS -q pub_gpu
#PBS -l walltime=400:00:00
#PBS -j oe

cd $PBS_O_WORKDIR
nvidia-smi -l 1000 > gpuStat.log 2>&1 &
source /public/spst/home/hupj/.bashrc
conda activate nequip

export WANDB_MODE=offline
export CUDA_VISIBLE_DEVICES=1
nequip-train nequipFull.yaml
