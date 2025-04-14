#!/bin/bash
#SBATCH -J vasp544  #作业名
#SBATCH -p ihicnormal #队列名
#SBATCH -N 1 #节点数
#SBATCH --ntasks-per-node=28     #节点进程数
#SBATCH --cpus-per-task=1     #每进程占用核心数
##SBATCH --exclusive  ##独占节点，按节点计费
#SBATCH -o %j.out
#SBATCH -e %j.err
                                                                                                                                                     

NPROCS=$SLURM_JOB_CPUS_PER_NODE

module purge; module load compiler/intel/2017.5.239  mpi/intelmpi/2017.4.239
export MKL_DEBUG_CPU_TYPE=5 #加速代码

export MKL_CBWR=AVX2 #使cpu默认支持avx2

export I_MPI_PIN_DOMAIN=numa #内存位置与cpu位置绑定，加速内存读取。对于内存带宽要求高的计算提速明显

export PATH=/public/home/ac877eihwp/software/vasp.5.4.4all/bin:$PATH #vasp程序路径

ulimit -s unlimited
for i in {16..20}; do
        echo -en "\033[38;5;${i}m###\033[0m"
done
echo -en " \033[38;5;87mAUTOEXEC SCRIPT RUNNING LOG\033[0m "
for i in {20..16}; do
        echo -en "\033[38;5;${i}m###\033[0m"
done; echo # skip to new line

if [[ ! -d "finished" ]]; then
        mkdir finished
fi
for i in `ls -d */`; do
        if [[ $i == "finished/" || ! -d $i ]]; then
                continue
        fi
        cd $i
        echo "`date +'[%Y/%m/%d %H:%M:%S]'` Starting running job in $i ..."
        srun --mpi=pmi2 -n $NPROCS vasp_std > resLog.out 2>&1
        echo "`date +'[%Y/%m/%d %H:%M:%S]'` Executed job in $i successfully"
        cd ..
        mv $i finished
done
echo "`date +'[%Y/%m/%d %H:%M:%S]'` All jobs done"
