# Workflow of building machine learning potentials via REVICO sampling

*REVICO (Random Exploration Via Imaginary Compounds Optimization) is a sampling method developed by Changxi Yang*

*Tutorial version 0.1 is written by Chenyu Wu, 2024/05/24*

---

**NOTE**: This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

---

## 0. Basic knowledge & environment configuration

Before starting, users should ensure that one can understand the following commands:

    $ conda env list
    $ conda activate <ENV_NAME>
    $ pip list

and know the catastrophic consequence of the following commands:

    $ rm -rf $tmpFolder/*       # var `tmpFolder` undefined
    $ cd $myPath && rm -rf *    # var `myPath` undefined
    $ sudo chmod 000 /*

Otherwise, please STOP USING THIS PACKAGE IMMEDIATELY.


The following programs & packages are needed for REVICO workflow:

&emsp;- python >= 3.10

&emsp;- nequip >= 0.5.6

&emsp;- numpy >= 1.25.2, numba >= 0.58.1, ase >= 3.22.1, scikit-learn >= 1.3.2, dscribe >= 2.1.0


## 1. Chemical elements selection & imaginary compounds generation

No prior knowledge except elements of chemical systems is needed.


#### Recommended count of imaginary compounds for several systems:

&emsp;- Au,N,C,H,O: 5000

<!-- &emsp;- Ag,Pd,C,H,O: 5000 -->

<!-- &emsp;- Pd,C,H: 3000 -->

<!-- &emsp;- Pt,O: 1000 -->


#### scripts:

&emsp;- randomGenerator.py *(not completed)*


#### command & usage:

&emsp;- Count structures in extxyz file:

        $ grep -c Latt STRUCTURE_FILES[.xyz]


## 2. Imaginary compounds relaxation via low precision DFT

Part of imaginary compounds generated in the previous step will be sampled for relaxation via low precision DFT.

More structures can be generated and more configuration space will be explored via low precision DFT with less computational cost.

Energies and forces obtained at this step should not be employed for final potential training, only for prototype model training.

Configurations with unconverged SCF calculation should be purged.


#### Recommended modification of DFT settings for low precision DFT:

&emsp;- ENCUT: 450 eV -> 300 eV

&emsp;- NSW: 0 -> 2000

&emsp;- KPOINTS: auto 25 -> Gamma

&emsp;- VASP version: vasp_std -> vasp_gam


#### Recommended count of sampled imaginary compounds for relaxation for several systems:

&emsp;- Au,N,C,H,O: 300


#### scripts:

&emsp;- INCAR_lowprec, KPOINTS_lowprec

&emsp;- batchVaspPrepare.py

&emsp;- jobSubmit_serial.sh

&emsp;- selectConvergedSCF.py


#### command & usage:

&emsp;- Prepare VASP files for parallel tasks:

        $ python batchVaspPrepare.py STRUCTURE_FILES[.xyz,...] BATCH_SIZE

&emsp;&emsp; For example, prepare parallel VASP tasks for random_PdCH.xyz in example/batchVaspPrepare:

        $ cd example/batchVaspPrepare
        $ grep -c Latt random_PdCH.xyz                                     # 100 structures in random_PdCH.xyz
        $ conda activate nequip                                            # python environment with ASE is needed
        $ cp ../../scripts/batchVaspPrepare.py ./
        $ cp ../../scripts/INCAR_lowprec INCAR
        $ cp ../../scripts/KPOINTS_lowprec KPOINTS
        $ cp /public/spst/home/hupj/apps/vasp/potpaw_PBE/Pd/POTCAR POT_Pd  # rename POTCAR of <ELEMENT> to POT_<ELEMENT>
        $ cp /public/spst/home/hupj/apps/vasp/potpaw_PBE/C/POTCAR POT_C
        $ cp /public/spst/home/hupj/apps/vasp/potpaw_PBE/H/POTCAR POT_H
        $ python batchVaspPrepare.py random_PdCH.xyz 50                    # split into 2 parallel tasks
        $ cp ../../scripts/jobSubmit_serial.sh 1_batch
        $ cp ../../scripts/jobSubmit_serial.sh 2_batch
        $ cd 1_batch; qsub jobSubmit_serial.sh; cd ..                      # submit 2 parallel tasks to PBS
        $ cd 2_batch; qsub jobSubmit_serial.sh; cd ..

&emsp;- Check state of SCF convergence manually by greping stdout of VASP calculation:

        $ grep elec resLog.out

&emsp;- Merge results of VASP calculations with unconverged SCF calculation purged:

        $ python selectConvergedSCF.py FOLDER_CONTAINS_OUTCAR

&emsp;&emsp; For example, merge results and purge unconverged SCF calculations in example/selectConvergedSCF/1_batch:

        $ conda activate nequip                         # python environment with ASE is needed
        $ cp ../../../scripts/selectConvergedSCF.py ./
        $ for i in `ls -d */`; do
        >     python selectConvergedSCF.py $i
        > done > resLog.out 2>&1                        # save running log to resLog.out
        $ grep -c Latt merged.xyz                       # check structures in merged.xyz


## 3. Data cleaning & selection for representative configurations

SST (Structures Select Twice) workflow can be employed for selection for representative configurations.

Dataset for validation and test will be generated out of selected configurations.


#### Recommended count of representative configurations for several systems:

&emsp;- 5 elements: 15,000


#### scripts:

&emsp;- SST_v0.2 package

&emsp;- jobSubmit_SST.sh


#### command & usage:

&emsp;- Launching SST workflow to obtain dataset for training and validation 
(parameter `elements`, `clsLimit`, `redCount` in the last two lines of main.py should be set before running):

        $ conda activate nequip
        $ python main.py source.xyz

&emsp;&emsp; **Tips**: SST tasks should be submitted to exec hosts if an out-of-memory error occurred. Use script jobSubmit_SST.sh.


## 4. Prototype model training

E(3)-equivariant graph neural network NequIP is employed to couple this workflow.

The prototype model will be used as optimizer of random configurations in large scale sampling.


#### scripts:

&emsp;- nequipFull.yaml

&emsp;- jobSubmit_nequip.sh


#### command & usage:

&emsp;- Necessary modification of nequipFull.yaml:

&emsp;&emsp; `run_name`, `dataset_file_name`, `chemical_symbols`, `validation_dataset_file_name`, `n_train`, `n_val` should be reset in a new training task.

&emsp;- Check GPU state while training:

&emsp;&emsp; Command `nvidia-smi` in jobSubmit_nequip.sh will record GPU state in file gpuStat.log while training.

&emsp;&emsp; One can use environment variable `CUDA_VISIBLE_DEVICES` to alter GPU for training on multi-card nodes.

&emsp;- Deploy a trained model:

        $ conda activate nequip
        $ nequip-deploy build --train-dir /path/to/checkpoints/of/model prototypeModel.pth

&emsp;&emsp; Here, file prototypeModel.pth is the deployed trained prototype model.


## 5. Large scale imaginary compounds relaxation via prototype machine learning potential

More structures can be generated and more configuration space will be explored via prototype machine learning potential with less computational cost.

Energies and forces obtained at this step should not be employed for final potential training, only for further data selection.

Configurations with abnormal energy or forces should be purged.


#### Recommended count of imaginary compounds for relaxation for several systems:

&emsp;- Au,N,C,H,O: 4700 (unused part of total 5000 configurations)


#### scripts:

&emsp;- nequipSampleViaOpt.py

&emsp;- jobSubmit_potSample.sh


#### command & usage:

&emsp;- 


## 6. Data cleaning & selection for representative configurations

Like step 3, SST workflow will be employed for data cleaning and selection.


#### Recommended count of representative configurations for several systems:

&emsp;- 5 elements: ~2% of total data (no less than 25,000)


#### scripts:

&emsp;- SST_v0.2 package

&emsp;- jobSubmit_SST.sh


#### command & usage:

&emsp;Same with step 3.


## 7. Labeling via high precision DFT

Energies and forces obtained at this step are the labels for final training.

Configurations with unconverged SCF calculation should be purged.


#### scripts:

&emsp;- INCAR_highprec, KPOINTS_highprec

&emsp;- batchVaspPrepare.py

&emsp;- jobSubmit_serial.sh

&emsp;- selectConvergedSCF.py


#### command & usage:

&emsp;Same with step 2.


## 8. Final model training

E(3)-equivariant graph neural network NequIP is employed to couple this workflow.

Final model will be obtained with configurations labeled via high precision DFT.


#### scripts:

&emsp;- nequipFull.yaml

&emsp;- jobSubmit_nequip.sh


#### command & usage:

&emsp;Same with step 4.

