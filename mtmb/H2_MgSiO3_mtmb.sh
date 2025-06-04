#based on a previous MLP (Peng & Deng, 2024 GRL). Our goal is to expand the explored H2 concentration.
#consider six H2 concentration, 32MgSiO3 + 16/32/48/64/80/96H2. 
#1-60 GPa, 2000-8000 K.


# (1) first run six scanned npt simulations in lammps using the six H2 concentration.
      # Design npt scan simulations. Different scan paths.

# (2) use asap to do PCA analysis for each dump file and extract 100 frames from each dump file.
asap gen_desc -s 1 --fxyz npt.dump soap -e -c 6 -n 4 -l 4 -g 0.44
python /home/hl0876/script/mldp/asap/select_frames.py -i ASAP-desc.xyz -n 100
mkdir pre
cd pre
python /home/hl0876/script/mldp/extract_deepmd.py -f ../npt.dump -fmt dump -id ../test-frame-select-fps-n-100.index -st -t 5000 #change the temperature 5000 K
#change type_map.raw


# (3) recalculate these frames using VASP, test MLP's performance, add some selected frames to train dataset, add some other frames to test dataset. 
mkdir inputs
cd inputs
#creat INCAR, KPOINTS, POTCAR, sub_vasp.sh
cd ..

#above code can be employed using npt_extract.sh

python /home/hl0876/script/mldp/recal_dpdata.py -d deepmd/ -if ./inputs/ -sc sbatch
cd recal 
#python /home/hl0876/script/mldp/post_recal_v2.py
#python /home/hl0876/script/mldp/check_nbands_nelm.py -ip all -v
# change job time if the calculation is not done, change NBAND value in the INCAR if bad band. If change the INCAR in the input directory, need to copy them to respective frame folder.
python /home/hl0876/script/mldp/post_recal_rerun.py -ip all -v -ss ../../inputs/sub_vasp.sh
source rerun
python /home/hl0876/script/mldp/merge_out.py -o OUTCAR -r y
rm -rfv ../deepmd
python /home/hl0876/script/mldp/extract_deepmd.py -d deepmd -ttr 10000
dp test -m /scratch/gpfs/JIEDENG/Haiyang/H2-MgSiO3/MTMB/Peng_Jie_raw/machine_learning_potential/m12v3.comp.pb -d m12v3
python /home/hl0876/script/mldp/model_dev/analysis.py -tf . -mp m12v3 -rf . -euc 0.12 -fuc 4 -flc 0.4 -elc 0.012
rm -rfv deepmd
# Initially use e_and_f so that the criterion is more strict and later (last a few iterations) use e_or_f so that you can select some frames
python /home/hl0876/script/mldp/extract_deepmd.py -f ./OUTCAR -id ./m12v3_id_e_and_f
























































































