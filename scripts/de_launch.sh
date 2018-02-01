#!/bin/bash

#PBS -N PP_RNAseq_DE_analysis
#PBS -k o
#PBS -l nodes=1:ppn=16,vmem=80gb
#PBS -l walltime=2:00:00
#PBS -m abe

module load r

WD=/N/dc2/scratch/rtraborn/T502_RNAseq

cd $WD

echo "Making symbolic links to bam files"

mkdir seudDir
cd seudDir

#ln -s ##### add symbolic links to bams

mkdir NHR40DIR
cd NHR40DIR
#ln -s ##### add symbolic links to bams

cd ../scripts

echo "Launching DE expression job"

R CMD BATCH de_analysis_neuron_mm10.R

echo "Done."

exit
