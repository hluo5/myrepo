#!/bin/bash

# Define the base path for the source and target locations
source_base="../../../it7_mtmb"
target_base="."

# List of the folder names
h2_folders=("32MgSiO3_16H2" "32MgSiO3_32H2" "32MgSiO3_48H2" "32MgSiO3_64H2" "32MgSiO3_80H2" "32MgSiO3_96H2")
#temperature_folders=("2000K" "3000K" "4000K" "5000K" "6000K" "7000K" "8000K")
temperature_folders=("1-30GPa_2000-5000K" "1-30GPa_5000-8000K" "30-60GPa_2000-5000K" "30-60GPa_5000-8000K")

# Loop through each folder combination
for h2_folder in "${h2_folders[@]}"; do
  for temp_folder in "${temperature_folders[@]}"; do
    # Construct the source file path
    source_file="$source_base/$h2_folder/$temp_folder/pre/recal/deepmd"
    
    # Construct the target file path
    target_file="$target_base/Luo.$h2_folder.$temp_folder.mtmb7.recal"
    
    # Copy the file
    cp -r "$source_file" "$target_file"
    
    echo "Copied from $source_file to $target_file"
  done
done

