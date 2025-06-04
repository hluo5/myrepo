#!/bin/bash
#SBATCH --job-name=simila         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1              # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:30:00          # total run time limit (HH:MM:SS)

module purge
module load anaconda3/2024.6
conda activate deepmd


mkdir similarity
cd similarity
#cp ../conf.lmp .
#cp ../../plumed.dat .
cp ../../../in.lammps .

lmp -in in.lammps
python /home/hl0876/script/mldp/similarity/merge_mgsiofe.py 

mkdir nw2
cd nw2
python ~/script/mldp/similarity/stat.py -n 2 -m mass -f ../merge.xyz

#cd ..

#mkdir nw0
#cd nw0
#python ~/script/mldp/similarity/stat.py -n 0 -m mass -f ../merge.xyz

#cd ..

#mkdir nw1
#cd nw1
#python ~/script/mldp/similarity/stat.py -n 1 -m mass -f ../merge.xyz 


