#!/bin/bash

set -e

WDIR=$(pwd)
MLDP=/home/hl0876/script/mldp

#for P in 50 250 500 1000 1500
#do
#for T in 2000 3000 4000 5000 6000 7000 8000
for T in 2000 3000 4000 5000
do

cd $WDIR/${T}K

cat > dft_select.sh << EOF
#!/bin/bash
#SBATCH --job-name=${T}K         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=12G        # memory per cpu-core (4G is default)
#SBATCH --time=00:05:00          # total run time limit (HH:MM:SS)

module purge
module load anaconda3/2024.6
conda activate deepmd

cd pre
cd recal
python $MLDP/merge_out.py -o OUTCAR -r y
rm -rfv ../deepmd
python $MLDP/extract_deepmd.py -d deepmd -ttr 10000
dp test -m ../../../../MLP_used/npt3_compressed.pb -d npt3
python $MLDP/model_dev/analysis.py -tf . -mp npt3 -rf . -euc 0.12 -fuc 4 -flc 0.4 -elc 0.012
rm -rfv deepmd
python $MLDP/extract_deepmd.py -f ./OUTCAR -id ./npt3_id_e_and_f
cd deepmd
mkdir training_data
mkdir validation_data
cp -r set.000 ./training_data
cp -r set.001 ./validation_data
cp type_map.raw ./training_data
cp type_map.raw ./validation_data
cp type.raw ./training_data
cp type.raw ./validation_data
rm -rfv set.000 set.001 type_map.raw type.raw



EOF

sbatch dft_select.sh

done
