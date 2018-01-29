#!/bin/bash

#PBS -N PP_RNAseq_fastqc_batch_T502
#PBS -k o
#PBS -l nodes=1:ppn=16,vmem=32gb
#PBS -l walltime=1:00:00
#PBS -m abe

module load fastqc

fileDir=/N/dc2/scratch/rtraborn/T502_fastqs/PP_RNAseq
####### Before running the script, please enter path to desired output directory, below ####
fqDir=<provide path to desired output directory>

cd $fqDir

ln -s ${fileDir}/GSF1659-NHR40-1_S5_R1_001.fastq.gz NHR40-1.R1.fastq.gz
ln -s ${fileDir}/GSF1659-NHR40-1_S5_R2_001.fastq.gz NHR40-1.R2.fastq.gz
ln -s ${fileDir}/GSF1659-NHR40-2_S6_R1_001.fastq.gz NHR40-2.R1.fastq.gz
ln -s ${fileDir}/GSF1659-NHR40-2_S6_R2_001.fastq.gz NHR40-2.R2.fastq.gz
ln -s ${fileDir}/GSF1659-NHR40-3_S7_R1_001.fastq.gz NHR40-3.R1.fastq.gz
ln -s ${fileDir}/GSF1659-NHR40-3_S7_R2_001.fastq.gz NHR40-3.R2.fastq.gz
ln -s ${fileDir}/GSF1659-NHR40-4_S8_R1_001.fastq.gz NHR40-4.R1.fastq.gz
ln -s ${fileDir}/GSF1659-NHR40-4_S8_R2_001.fastq.gz NHR40-4.R2.fastq.gz

ln -s ${fileDir}/GSF1659-Seud1-1_S1_R1_001.fastq.gz Seud1-1.R1.fastq.gz
ln -s ${fileDir}/GSF1659-Seud1-1_S1_R2_001.fastq.gz Seud1-1.R2.fastq.gz
ln -s ${fileDir}/GSF1659-Seud1-2_S2_R1_001.fastq.gz Seud1-2.R1.fastq.gz
ln -s ${fileDir}/GSF1659-Seud1-2_S2_R2_001.fastq.gz Sedu1-2.R2.fastq.gz
ln -s ${fileDir}/GSF1659-Seud1-3_S3_R1_001.fastq.gz Seud1-3.R1.fastq.gz
ln -s ${fileDir}/GSF1659-Seud1-3_S3_R2_001.fastq.gz Seud1-3.R2.fastq.gz
ln -s ${fileDir}/GSF1659-Seud1-4_S4_R1_001.fastq.gz Seud1-4.R1.fastq.gz
ln -s ${fileDir}/GSF1659-Seud1-4_S4_R2_001.fastq.gz Seud1-4.R2.fastq.gz

echo "Starting job"

echo "Running fastqc on the RNA-seq files"

for fq in *.fastq.gz; do
echo "Starting fastqc for $fq"
fastqc $fq
done

echo "Fastqc job is complete"

exit


