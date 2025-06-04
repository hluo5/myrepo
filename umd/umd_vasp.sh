#!/bin/bash
#SBATCH --job-name=umd        # create a short name for your job
#SBATCH --nodes=1                 # node count
#SBATCH --ntasks-per-node=1      # total number of tasks per node
#SBATCH --cpus-per-task=1         # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=2G         # memory per cpu-core (4G is default)
#SBATCH --time=00:10:00           # total run time limit (HH:MM:SS)
#SBATCH --account=jiedeng

#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
#export PATH=$PATH:/scratch/gpfs/JIEDENG/Haiyang/software/vasp.6.4.2/bin

#UMDcommand
module purge
module load anaconda3/2024.6
conda activate myenv

#cp XDATCAR XDAT
#tail -n +8 XDATCAR > XDATCAR0
#sed '/Direct/d' ./XDATCAR0 > XDATCAR1
#rm XDATCAR0
#mv XDATCAR1 XDATCAR



#New UMD version

python /scratch/gpfs/JIEDENG/Haiyang/software/MAGMATOMIX-UMD_package-2.2/src/VaspParser.py -f OUTCAR -i 5000 >> UMD.outlog

#bash /home/hl0876/bin/block >> UMD.outlog

python /scratch/gpfs/JIEDENG/Haiyang/software/MAGMATOMIX-UMD_package-2.2/src/fullaverages.py -u 0 >> UMD.outlog

python /scratch/gpfs/JIEDENG/Haiyang/software/MAGMATOMIX-UMD_package-2.2/src/gofr_umd.py -f OUTCAR.umd.dat -s 1 -d 0.01 >> UMD.outlog

#python /scratch/gpfs/JIEDENG/Haiyang/software/MAGMATOMIX-UMD_package-2.2/src/bonding_umd.py -f OUTCAR.umd.dat -s 1 -i bonds.input >> UMD.outlog

#python /scratch/gpfs/JIEDENG/Haiyang/software/MAGMATOMIX-UMD_package-2.2/src/speciation_and_angles.py -b OUTCAR.bonding.dat -f OUTCAR.umd.dat -s 1 -c Si -a O -r 1 -t 1 -l 1 >> UMD.outlog

#python /scratch/gpfs/JIEDENG/Haiyang/software/MAGMATOMIX-UMD_package-2.2/src/msd_umd.py -f OUTCAR.umd.dat -z 50 -v 50 >> UMD.outlog

#python /scratch/gpfs/JIEDENG/Haiyang/software/MAGMATOMIX-UMD_package-2.2/src/vibr_spectrum_umd_fast.py -f OUTCAR.umd.dat -t 3500 >> UMD.outlog

#python /scratch/gpfs/JIEDENG/Haiyang/software/MAGMATOMIX-UMD_package-2.2/src/viscosity_umd.py -f OUTCAR.umd.dat >> UMD.outlog

