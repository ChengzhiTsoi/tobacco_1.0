#!/bin/bash

#$ -N eqeq
#$ -q all.q@apollo3-c19
#$ -cwd

topdir=$PWD
N_samp=$(wc -l < MOFs_gcmc.csv) 

sed -i 's/\r$//' MOFs_gcmc.csv
for ((i=1; i<N_samp; i++))
do
    Row=$((i + 1))
    read FrameworkName <<< $(awk 'FNR=='$Row' {print $1}' MOFs_gcmc.csv)

    cat ${topdir}/EQeq/FrameworkName_EQeq_Ewald_1.20_-2.00.mol >> ${topdir}/prepared_mofs/FrameworkName.mol   

done
