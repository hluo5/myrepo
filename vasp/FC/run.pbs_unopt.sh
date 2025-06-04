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
10.436951000000000      0.0000000000000000      0.0000000000000000 
 0.000000000000000     10.4369510000000000      0.0000000000000000
 0.000000000000000      0.0000000000000000     10.4369510000000000
Mg Si O
16 16 48
`sed -n "$(grep -n "Direct configuration=  $i" XDAT | tail -1 | cut -d : -f 1),+80p" XDAT`
!
`sed -i '7a Selective Dynamics' POSCAR`
`sed -i '10s/$/ T T T/' POSCAR`
`sed -i '11,89s/$/ F F F/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 1 $FX $Fy $FZ >>CFC.txt

`sed -i '10s/T T T/F F F/' POSCAR`
`sed -i '11s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 2 $FX $Fy $FZ >>CFC.txt

`sed -i '11s/T T T/F F F/' POSCAR`
`sed -i '12s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 3 $FX $Fy $FZ >>CFC.txt

`sed -i '12s/T T T/F F F/' POSCAR`
`sed -i '13s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 4 $FX $Fy $FZ >>CFC.txt

`sed -i '13s/T T T/F F F/' POSCAR`
`sed -i '14s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 5 $FX $Fy $FZ >>CFC.txt

`sed -i '14s/T T T/F F F/' POSCAR`
`sed -i '15s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 6  $FX $Fy $FZ >>CFC.txt

`sed -i '15s/T T T/F F F/' POSCAR`
`sed -i '16s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 7 $FX $Fy $FZ >>CFC.txt

`sed -i '16s/T T T/F F F/' POSCAR`
`sed -i '17s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 8 $FX $Fy $FZ >>CFC.txt

`sed -i '17s/T T T/F F F/' POSCAR`
`sed -i '18s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 9 $FX $Fy $FZ >>CFC.txt

`sed -i '18s/T T T/F F F/' POSCAR`
`sed -i '19s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 10 $FX $Fy $FZ >>CFC.txt

`sed -i '19s/T T T/F F F/' POSCAR`
`sed -i '20s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 11 $FX $Fy $FZ >>CFC.txt

`sed -i '20s/T T T/F F F/' POSCAR`
`sed -i '21s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 12 $FX $Fy $FZ >>CFC.txt

`sed -i '21s/T T T/F F F/' POSCAR`
`sed -i '22s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 13 $FX $Fy $FZ >>CFC.txt

`sed -i '22s/T T T/F F F/' POSCAR`
`sed -i '23s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 14 $FX $Fy $FZ >>CFC.txt

`sed -i '23s/T T T/F F F/' POSCAR`
`sed -i '24s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 15 $FX $Fy $FZ >>CFC.txt

`sed -i '24s/T T T/F F F/' POSCAR`
`sed -i '25s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 16 $FX $Fy $FZ >>CFC.txt

`sed -i '25s/T T T/F F F/' POSCAR`
`sed -i '26s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 17 $FX $Fy $FZ >>CFC.txt

`sed -i '26s/T T T/F F F/' POSCAR`
`sed -i '27s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 18 $FX $Fy $FZ >>CFC.txt

`sed -i '27s/T T T/F F F/' POSCAR`
`sed -i '28s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 19 $FX $Fy $FZ >>CFC.txt

`sed -i '28s/T T T/F F F/' POSCAR`
`sed -i '29s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 20 $FX $Fy $FZ >>CFC.txt

`sed -i '29s/T T T/F F F/' POSCAR`
`sed -i '30s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 21 $FX $Fy $FZ >>CFC.txt

`sed -i '30s/T T T/F F F/' POSCAR`
`sed -i '31s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 22 $FX $Fy $FZ >>CFC.txt

`sed -i '31s/T T T/F F F/' POSCAR`
`sed -i '32s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 23 $FX $Fy $FZ >>CFC.txt

`sed -i '32s/T T T/F F F/' POSCAR`
`sed -i '33s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 24 $FX $Fy $FZ >>CFC.txt

`sed -i '33s/T T T/F F F/' POSCAR`
`sed -i '34s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 25 $FX $Fy $FZ >>CFC.txt

`sed -i '34s/T T T/F F F/' POSCAR`
`sed -i '35s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 26 $FX $Fy $FZ >>CFC.txt

`sed -i '35s/T T T/F F F/' POSCAR`
`sed -i '36s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 27 $FX $Fy $FZ >>CFC.txt

`sed -i '36s/T T T/F F F/' POSCAR`
`sed -i '37s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 28 $FX $Fy $FZ >>CFC.txt

`sed -i '37s/T T T/F F F/' POSCAR`
`sed -i '38s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 29 $FX $Fy $FZ >>CFC.txt

`sed -i '38s/T T T/F F F/' POSCAR`
`sed -i '39s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 30 $FX $Fy $FZ >>CFC.txt

`sed -i '39s/T T T/F F F/' POSCAR`
`sed -i '40s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 31 $FX $Fy $FZ >>CFC.txt

`sed -i '40s/T T T/F F F/' POSCAR`
`sed -i '41s/F F F/T T T/' POSCAR`
yhrun -N $SLURM_NNODES -n $SLURM_NTASKS ${VM} >& log.$SLURM_JOBID
FX=`grep "SECOND" -A3 OUTCAR | tail -1 | cut -c 7-18`
Fy=`grep "SECOND" -A4 OUTCAR | tail -1 | cut -c 19-30`
FZ=`grep "SECOND" -A5 OUTCAR | tail -1 | cut -c 31-41`
echo $i 32 $FX $Fy $FZ >>CFC.txt
done









