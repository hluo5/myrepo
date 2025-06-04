#!/bin/bash
  
#SBATCH --job-name=index800
#SBATCH --nodes=1                #node count
#SBATCH --ntasks-per-node=1  #total number of tasks per node
#SBATCH --cpus-per-task=1  #cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=2G
#SBATCH --time=1:00:00      
#SBATCH --constraint=cascade

module purge
module load anaconda3/2021.11
conda activate dpdev

mkdir index_800
cd index_800
asap gen_desc -s 1 --fxyz ../OUTCAR soap -e -c 6 -n 6 -l 6 -g 0.44
#asap gen_desc -s 2 --fxyz npt.dump soap -e -c 6 -n 4 -l 4 -g 0.44
python /home/hl0876/script/mldp/asap/select_frames.py -i ASAP-desc.xyz -n 800
mkdir pre
cd pre
python /home/hl0876/script/mldp/extract_deepmd.py -f ../../OUTCAR -id ../test-frame-select-fps-n-800.index -st # OUTCAR contains temperature info
#python /home/hl0876/script/mldp/extract_deepmd.py -f ../npt.dump -fmt dump -id ../test-frame-select-fps-n-100.index -st -t 4000

#prepare a folder named inputs including INCAR KPOINTS POTCAR sub_vasp.sh

#python /home/hl0876/script/mldp/recal_dpdata.py -d deepmd/ -if ../../inputs/ -sc sbatch

#python /home/hl0876/script/mldp/recal_dpdata.py -d deepmd/ -if ./inputs/ -rv no

#cd recal
#python /home/hl0876/script/mldp/check_nbands_nelm.py -ip all -v
#python /home/hl0876/script/mldp/merge_out.py -o OUTCAR -r y
#python /home/hl0876/script/mldp/extract_deepmd.py -d deepmd -ttr 10000

#start training mpv0 during the 0th round
#since the 1st round, do dp test and ......



