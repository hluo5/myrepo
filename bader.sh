#!/bin/bash

#SBATCH --job-name=1000K
#SBATCH --nodes=8                #node count
#SBATCH --ntasks-per-node=4  #total number of tasks per node
#SBATCH --cpus-per-task=1  #cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         #memory per cpu-core (4G per cpu-core is default)
#SBATCH --time=24:00:00      

#VASP6_tiger
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export PATH=$PATH:/home/hl0876/vasp.6.3.2/bin
#module purge
#module load intel/19.1/64/19.1.1.217 intel-mpi/intel/2019.1/64
#module load intel/19.0/64/19.0.5.281 intel-mpi/intel/2019.5/64
#srun --distribution=block:block --hint=nomultithread vasp_gam

#VASP5_tiger
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PATH=$PATH:/home/hl0876/vasp.5.4.4.pl2/bin
module purge
module load intel/19.1/64/19.1.1.217 intel-mpi/intel/2019.7/64
#module load intel/19.0/64/19.0.5.281 intel-mpi/intel/2019.5/64
#srun --distribution=block:block --hint=nomultithread vasp_gam

export qvasppath=/home/hl0876/qvasp-v2.22
export PATH=/home/hl0876/qvasp-v2.22:$PATH
export PATH=/home/hl0876/qvasp-v2.22/exefile/Tools/USERTooLs/vtstscripts:$PATH

for i in $(seq 2000 250 6000)
do
cat >POSCAR <<!
MgFeSiO3
           1
     8.536860    0.000000    0.000000
     0.000000    8.536860    0.000000
     0.000000    0.000000    8.536860
   Mg   Fe   Si   O    K
  28   4  32  98   1
`sed -n "$(grep -n "Direct configuration=  $i$" XDAT | tail -1 | cut -d : -f 1),+163p" XDAT`
!
srun --distribution=block:block --hint=nomultithread vasp_gam
qvasp -baderd
#XVal=`sed -n '156p' ACF.dat | cut -c 44-54`
#XVol=`sed -n '156p' ACF.dat | cut -c 70-79`
#echo $i $XVal $XVol >> Val.txt
#done
PMVal=`awk 'NR==3,NR==165 {print $5}' ACF.dat`
PMVol=`awk 'NR==3,NR==165 {print $7}' ACF.dat`
echo $i $PMVal $PMVol >> PMVal.txt
done
awk '{s+=$2; ss+=$2^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' PMVal.txt >> PMVal_avg.txt
awk '{s+=$3; ss+=$3^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' PMVal.txt >> PMVal_avg.txt

