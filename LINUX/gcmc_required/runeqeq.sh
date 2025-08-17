o#!/bin/sh
topdir=$PWD
N_samp=$(wc -l < MOFs_gcmc.csv) 

sed -i 's/\r$//' MOFs_gcmc.csv
for ((i=1; i<N_samp; i++))
do
    Row=$((i + 1))
    read FrameworkName <<< $(awk 'FNR=='$Row' {print $1}' MOFs_gcmc.csv)

    cd ${topdir}/EQeq
    ./eqeq FrameworkName.cif    

    echo $COUNTER

done
