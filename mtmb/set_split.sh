#!/bin/bash

# Define the parent directory 'a'
parent_dir="train_test_data"

# Loop through each folder (A, B, C, D, etc.) in the parent directory
for folder in "$parent_dir"/*; do
  if [ -d "$folder" ]; then
    # Create training_data and validation_data directories if they don't exist
    mkdir -p "$folder/training_data"
    mkdir -p "$folder/validation_data"
    
    # Find the folder that starts with 'deepmd' (and may have extra characters)
    deepmd_dir=$(find "$folder" -type d -name "deepmd*" | head -n 1)
    
    if [ -d "$deepmd_dir" ]; then
      # Copy the files for training_data (set.000, type_map.raw, type.raw)
      cp -r "$deepmd_dir/set.000" "$folder/training_data/"
      cp "$deepmd_dir/type_map.raw" "$folder/training_data/"
      cp "$deepmd_dir/type.raw" "$folder/training_data/"
      
      # Copy the files for validation_data (set.001, type_map.raw, type.raw)
      cp -r "$deepmd_dir/set.001" "$folder/validation_data/"
      cp "$deepmd_dir/type_map.raw" "$folder/validation_data/"
      cp "$deepmd_dir/type.raw" "$folder/validation_data/"
      
      # Delete the 'deepmd' folder (or any folder that starts with 'deepmd')
      rm -rf "$deepmd_dir"
    fi
  fi
done

