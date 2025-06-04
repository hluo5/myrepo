#!/bin/bash

rm PMVal_avg.txt
for i in $(seq 2 1 29)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> Mg_Val.txt
done
for i in $(seq 30 1 33)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> Fe_Val.txt
done
for i in $(seq 34 1 65)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> Si_Val.txt
done
for i in $(seq 66 1 163)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> O_Val.txt
done
for i in $(seq 164 1 164)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> K_Val.txt
done
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' Mg_Val.txt >> Mg_Val_avg.txt
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' Fe_Val.txt >> Fe_Val_avg.txt
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' Si_Val.txt >> Si_Val_avg.txt
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' O_Val.txt >> O_Val_avg.txt
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' K_Val.txt >> K_Val_avg.txt
rm Mg_Val.txt Fe_Val.txt Si_Val.txt O_Val.txt K_Val.txt

for i in $(seq 165 1 192)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> Mg_Vol.txt
done
for i in $(seq 193 1 196)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> Fe_Vol.txt
done
for i in $(seq 197 1 228)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> Si_Vol.txt
done
for i in $(seq 229 1 326)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> O_Vol.txt
done
for i in $(seq 327 1 327)
do
awk -v Val="$i" '{print $Val}' PMVal.txt >> K_Vol.txt
done                                                                       
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' Mg_Vol.txt >> Mg_Vol_avg.txt
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' Fe_Vol.txt >> Fe_Vol_avg.txt
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' Si_Vol.txt >> Si_Vol_avg.txt
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' O_Vol.txt >> O_Vol_avg.txt
awk '{s+=$1; ss+=$1^2} END{print m=s/NR, sqrt(ss/NR-m^2)/sqrt(NR)}' K_Vol.txt >> K_Vol_avg.txt
rm Mg_Vol.txt Fe_Vol.txt Si_Vol.txt O_Vol.txt K_Vol.txt
