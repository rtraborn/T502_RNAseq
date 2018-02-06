#!/bin/bash

#PBS -N PP_RNAseq_DE_analysis
#PBS -k o
#PBS -l nodes=1:ppn=16,vmem=80gb
#PBS -l walltime=2:00:00
#PBS -m abe

module load r

/<please enter the full path to scripts/>
#WD=/N/dc2/scratch/rtraborn/T502_RNAseq/scripts

cd $WD

echo "Launching DE expression job"

R CMD BATCH de_analysis_Prist_hybrid2.R

echo "Done."

exit
