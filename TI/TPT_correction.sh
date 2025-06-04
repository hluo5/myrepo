#!/bin/bash

cd ../0.0/01

T=3500

#vcaenergy_la1.0
#pressure5
grep 'F' OSZICAR | cut -c 37-51 > e.txt
grep 'total pressure' OUTCAR | cut -c 22-30 > P.txt

sed -n '1000~200p' e.txt > F_200
paste --delimiters=" " F_200  G600toK800/F.txt|awk '{if(NF>1) print $2-$1}' > F_diff
F_diff_avg=`awk '{s+=$1} END{print m=s/NR}' F_diff`
paste --delimiters=" " F_200  G600toK800/F.txt|awk -v F_diff_avg=$F_diff_avg '{if(NF>1) print ($2-$1-F_diff_avg)**2}' > F_diff2
F_diff_avg2=`awk '{s+=$1} END{print m=s/NR}' F_diff2`
KT=`echo $T | awk '{print $1*0.0000861733}'`
echo $F_diff_avg $F_diff_avg2 $KT | awk '{print $1-$2/2/$3}' > F_corr

sed -n '1000~200p' P_external.txt > P_200
paste --delimiters=" " P_200  G600toK800/P.txt|awk '{if(NF>1) print $2-$1}' > P_diff
awk '{s+=$1} END{print m=s/NR}' P_diff > P_diff_avg

cat F_corr
cat P_diff_avg
rm F_200 P_200 F_diff F_diff2 P_diff

cd ../../1.0/02

T=3500

#vcaenergy_la1.0
#pressure5
grep 'F' OSZICAR | cut -c 37-51 > e.txt
grep 'total pressure' OUTCAR | cut -c 22-30 > P.txt

sed -n '1000~200p' e.txt > F_200
paste --delimiters=" " F_200  G600toK800/F.txt|awk '{if(NF>1) print $2-$1}' > F_diff
F_diff_avg=`awk '{s+=$1} END{print m=s/NR}' F_diff`
paste --delimiters=" " F_200  G600toK800/F.txt|awk -v F_diff_avg=$F_diff_avg '{if(NF>1) print ($2-$1-F_diff_avg)**2}' > F_diff2
F_diff_avg2=`awk '{s+=$1} END{print m=s/NR}' F_diff2`
KT=`echo $T | awk '{print $1*0.0000861733}'`
echo $F_diff_avg $F_diff_avg2 $KT | awk '{print $1-$2/2/$3}' > F_corr

sed -n '1000~200p' P_external.txt > P_200
paste --delimiters=" " P_200  G600toK800/P.txt|awk '{if(NF>1) print $2-$1}' > P_diff
awk '{s+=$1} END{print m=s/NR}' P_diff > P_diff_avg

paste --delimiters=" " ../../0.0/01/F_corr  F_corr|awk '{if(NF>1) print $2-$1}' > F_corr_diff

cat F_corr
cat P_diff_avg
cat F_corr_diff
rm F_200 P_200 F_diff F_diff2 P_diff
