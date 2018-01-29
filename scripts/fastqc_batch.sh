#!/bin/bash

#PBS -N PP_RNAseq_fastqc_batch_T502
#PBS -k o
#PBS -l nodes=1:ppn=16,vmem=32gb
#PBS -l walltime=2:00:00
#PBS -m abe

module load fastqc

WD=/N/dc2/scratch/rtraborn/T502_fastqs/PP_RNAseq

cd $WD

echo "Starting job"

echo "Running fastqc on the RNA-seq files"

for fq in *.fastq.gz; do
echo "Starting fastqc for $fq"
fastqc $fq
done

echo "Fastqc job is complete"

exit


