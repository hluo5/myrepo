#!/bin/bash

rm PMVal_avg.txt
for i in $(seq 2 1 97)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> Mg_Val.txt
done
for i in $(seq 98 1 185)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> Si_Val.txt
done
for i in $(seq 186 1 457)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> O_Val.txt
done
for i in $(seq 458 1 545 )
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> Fe_Val.txt
done
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' Mg_Val.txt >> Mg_Val_avg.txt
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' Si_Val.txt >> Si_Val_avg.txt
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' O_Val.txt >> O_Val_avg.txt
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' Fe_Val.txt >> Fe_Val_avg.txt
#rm Mg_Val.txt Si_Val.txt O_Val.txt Fe_Val.txt 

for i in $(seq 546 1 641)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> Mg_Vol.txt
done
for i in $(seq 642 1 729)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> Si_Vol.txt
done
for i in $(seq 730 1 1001)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> O_Vol.txt
done
for i in $(seq 1002 1 1089)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> Fe_Vol.txt
done
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' Mg_Vol.txt >> Mg_Vol_avg.txt
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' Si_Vol.txt >> Si_Vol_avg.txt
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' O_Vol.txt >> O_Vol_avg.txt
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' Fe_Vol.txt >> Fe_Vol_avg.txt
#rm Mg_Vol.txt Si_Vol.txt O_Vol.txt Fe_Vol.txt 
