#!/bin/bash
  
#SBATCH --job-name=index
#SBATCH --nodes=1                #node count
#SBATCH --ntasks-per-node=1  #total number of tasks per node
#SBATCH --cpus-per-task=1  #cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00      

module purge
module load anaconda3/2024.6
conda activate deepmd

mkdir index_500
cd index_500
asap gen_desc -s 1 --fxyz ../OUTCAR soap -e -c 6 -n 6 -l 6 -g 0.44
python /home/hl0876/script/mldp/asap/select_frames.py -i ASAP-desc.xyz -n 500
python /home/hl0876/script/mldp/asap/select_frames.py -i ASAP-desc.xyz -n 1000
python /home/hl0876/script/mldp/asap/select_frames.py -i ASAP-desc.xyz -n 3000
python /home/hl0876/script/mldp/asap/select_frames.py -i ASAP-desc.xyz -n 5000
python /home/hl0876/script/mldp/asap/select_frames.py -i ASAP-desc.xyz -n 7000
python /home/hl0876/script/mldp/extract_deepmd.py -f ../OUTCAR -id ./test-frame-select-fps-n-500.index 
cd ..


mkdir index_1000
cd index_1000
cp ../index_500/ASAP-desc-state.yaml .
cp ../index_500/ASAP-desc.xyz .
cp ../index_500/test-frame-select-fps-n-1000.index .
python /home/hl0876/script/mldp/extract_deepmd.py -f ../OUTCAR -id ./test-frame-select-fps-n-1000.index 
cd ..

mkdir index_3000
cd index_3000
cp ../index_500/ASAP-desc-state.yaml .
cp ../index_500/ASAP-desc.xyz .
cp ../index_500/test-frame-select-fps-n-3000.index .
python /home/hl0876/script/mldp/extract_deepmd.py -f ../OUTCAR -id ./test-frame-select-fps-n-3000.index 
cd ..

mkdir index_5000
cd index_5000
cp ../index_500/ASAP-desc-state.yaml .
cp ../index_500/ASAP-desc.xyz .
cp ../index_500/test-frame-select-fps-n-5000.index .
python /home/hl0876/script/mldp/extract_deepmd.py -f ../OUTCAR -id ./test-frame-select-fps-n-5000.index
cd ..

mkdir index_7000
cd index_7000
cp ../index_500/ASAP-desc-state.yaml .
cp ../index_500/ASAP-desc.xyz .
cp ../index_500/test-frame-select-fps-n-7000.index .
python /home/hl0876/script/mldp/extract_deepmd.py -f ../OUTCAR -id ./test-frame-select-fps-n-7000.index 
cd ..


