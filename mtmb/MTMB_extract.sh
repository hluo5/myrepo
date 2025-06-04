#!/bin/bash

set -e

WDIR=$(pwd)
MLDP=/home/hl0876/script/mldp

# Group 1: Folders with temperature range 2000-5000K and T=3500
for FOLDER in 1-30GPa_2000-5000K 30-60GPa_2000-5000K; do
  T=3500
  cd "$WDIR/$FOLDER"
  fparam_value=$(awk "BEGIN {printf \"%.8f\", 8.617333262e-5 * $T}")

  cat > MTMB_extract.sh << EOF
#!/bin/bash
#SBATCH --job-name=${T}K         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)

module purge
module load anaconda3/2024.6
conda activate deepmd

asap gen_desc -s 1 --fxyz MTMB.dump soap -e -c 6 -n 4 -l 4 -g 0.44
python $MLDP/asap/select_frames.py -i ASAP-desc.xyz -n 50

mkdir pre
cd pre
python $MLDP/extract_deepmd.py -f ../MTMB.dump -fmt dump -id ../test-frame-select-fps-n-50.index -st -t $T

cd deepmd
sed -i "s/TYPE_0/Mg/" type_map.raw
sed -i "s/TYPE_1/Si/" type_map.raw
sed -i "s/TYPE_2/O/" type_map.raw
sed -i "s/TYPE_3/H/" type_map.raw
cd ..

mkdir inputs
cd inputs
cp ../../../inputs/INCAR .
cp ../../../inputs/KPOINTS .
cp ../../../inputs/POTCAR .
cp ../../../inputs/sub_vasp.sh .
sed -i "s/sigma/$fparam_value/" INCAR

EOF

  sbatch MTMB_extract.sh
done

# Group 2: Folders with temperature range 5000-8000K and T=6500
for FOLDER in 1-30GPa_5000-8000K 30-60GPa_5000-8000K; do
  T=6500
  cd "$WDIR/$FOLDER"
  fparam_value=$(awk "BEGIN {printf \"%.8f\", 8.617333262e-5 * $T}")

  cat > MTMB_extract.sh << EOF
#!/bin/bash
#SBATCH --job-name=${T}K         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)

module purge
module load anaconda3/2024.6
conda activate deepmd

asap gen_desc -s 1 --fxyz MTMB.dump soap -e -c 6 -n 4 -l 4 -g 0.44
python $MLDP/asap/select_frames.py -i ASAP-desc.xyz -n 50

mkdir pre
cd pre
python $MLDP/extract_deepmd.py -f ../MTMB.dump -fmt dump -id ../test-frame-select-fps-n-50.index -st -t $T

cd deepmd
sed -i "s/TYPE_0/Mg/" type_map.raw
sed -i "s/TYPE_1/Si/" type_map.raw
sed -i "s/TYPE_2/O/" type_map.raw
sed -i "s/TYPE_3/H/" type_map.raw
cd ..

mkdir inputs
cd inputs
cp ../../../inputs/INCAR .
cp ../../../inputs/KPOINTS .
cp ../../../inputs/POTCAR .
cp ../../../inputs/sub_vasp.sh .
sed -i "s/sigma/$fparam_value/" INCAR

EOF

  sbatch MTMB_extract.sh
done

