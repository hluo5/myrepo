#!/bin/bash
for i in 1 2 3
do
cd la$i
cd 01
grep 'free  energy   TOTEN  =' OUTCAR | awk '{print $5}' > F
cd ../02
grep 'free  energy   TOTEN  =' OUTCAR | awk '{print $5}' > F
cd ../
paste --delimiters=" "  01/F 02/F|awk '{if(NF>1) print $2-$1}' > F_diff
alchemyblock
cd ../
done
cat la1/F_diff_block.csv la2/F_diff_block.csv la3/F_diff_block.csv > F_diff_block.csv



