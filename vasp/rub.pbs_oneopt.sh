#!/bin/bash

######### VASP PATH #################
VASP_PATH=/HOME/ac_gig_hyluo/vasp5.3.3/bin_beta
#####################################

######### VASP MODE ##################
##
VM=vasp
######################################
export PATH=${VASP_PATH}:$PATH
export LD_LIBRARY_PATH=/HOME/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64:$LD_LIBRARY_PATH

for i in $(seq 5000 200 9800)
do
cat >POSCAR <<!
MgSiO3
1.000000000000
10.110580000000000      0.0000000000000000      0.0000000000000000 
 0.000000000000000     10.1105800000000000      0.0000000000000000
 0.000000000000000      0.0000000000000000     10.1105800000000000
Mg Si O
16 16 48
`sed -n "$(grep -n "Direct configuration=  $i" XDAT | tail -1 | cut -d : -f 1),+80p" XDAT`
!
`sed -i '7a Selective Dynamics' POSCAR`
`sed -i '10s/$/ T T T/' POSCAR`
`sed -i '11,89s/$/ F F F/' POSCAR`
cp POSCAR poscar
cp INCAR.youhua INCAR
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
cp CONTCAR POSCAR
cp INCAR.pinlv INCAR
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 1 $FX $Fy $FZ >>CFC.txt
 for j in $(seq 10 1 40)
 do
 k=$((j+1))
 m=$((j-8))
 `sed -i ''"$j"'s/T T T/F F F/' poscar`
 `sed -i ''"$k"'s/F F F/T T T/' poscar`
 cp poscar POSCAR
 cp INCAR.youhua INCAR
 yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
 cp CONTCAR POSCAR
 cp INCAR.pinlv INCAR
 yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
 FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
 Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
 FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
 echo $i $m $FX $Fy $FZ >>CFC.txt
 done
done










