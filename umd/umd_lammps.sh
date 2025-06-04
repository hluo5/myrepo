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

module purge
module load anaconda3/2024.6
conda activate myenv

#New UMD version
#to use LAMMPSParser2umd.py from non-zero steps, change the Reading data file ... (below serveral lines) and use deepmd-kit at:  

export UMD_PATH="/scratch/gpfs/JIEDENG/Haiyang/software/MAGMATOMIX-UMD_package/src"

python $UMD_PATH/modified_LAMMPSParser2umd.py -f npt.dump -l log.lammps -s 450 >> UMD.outlog

#bash /home/hl0876/bin/block >> UMD.outlog

#python $UMD_PATH/fullaverages.py -u 0 >> UMD.outlog

python $UMD_PATH/gofr_umd.py -f npt.dump.umd.dat -s 1 -d 0.01 >> UMD.outlog
#if jump frames using -s 10, pay attention to whether times the g(r) results by 10

#python $UMD_PATH/bonding_umd.py -f npt.dump.umd.dat -s 1 -i bonds.input >> UMD.outlog

#python $UMD_PATH/speciation_and_angles.py -b npt.dump.bonding.dat -f npt.dump.umd.dat -s 1 -c Si -a O -r 1 -t 1 -l 1 >> UMD.outlog

#python $UMD_PATH/msd_umd.py -f npt.dump.umd.dat -z 50 -v 50 >> UMD.outlog

#python $UMD_PATH/vibr_spectrum_umd_fast.py -f npt.dump.umd.dat -t 3500 >> UMD.outlog

#python $UMD_PATH/viscosity_umd.py -f npt.dump.umd.dat >> UMD.outlog
