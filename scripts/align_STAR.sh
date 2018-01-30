#!/bin/bash

#PBS -N PP_RNAseq_STAR_align_T502
#PBS -k o
#PBS -q debug
#PBS -l nodes=1:ppn=16,vmem=80gb
#PBS -l walltime=1:00:00
#PBS -m abe

module load star

fileDir=/N/dc2/scratch/rtraborn/T502_fastqs/PP_RNAseq
####### Before running the script, please enter path to desired output directory, below ####
WD=/N/u/rtraborn/Carbonate/T502_RNAseq
#WD=/N/u/<yourUserId>/Carbonate/T502_RNAseq/
outDir=/N/dc2/scratch/rtraborn/starAlign
fqDir=fastqs
genomedir=${WD}/fasta
genomeFasta=pacificus_Hybrid2.fa
configfile=${WD}/scripts/STARalign.conf
nThreads=16


function readconfigfile {
# Read the specified ($1) STARalign configuration file:
if [ ! -e "$1" ] ; then
  echo ""
  echo "Fatal error: STARalign config file $1 does not exist. Please check."
  exit 1
fi

STARgenomeGenerateOptions=`grep '^STARgenomeGenerateOptions=' "$1" | awk -F"=" '{print $2}'`
STARalignReadsOptions=`grep '^STARalignReadsOptions=' "$1" | awk -F"=" '{print $2}'`
bamFilterOptions=`grep '^bamFilterOptions=' "$1" | awk -F"=" '{print $2}'`
}

echo "Starting job"

echo "Making symbolic links to fastq files."

cd $WD

if [ ! -d "$fqDir" ] ; then
    mkdir $fqDir
  fi

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

cd $WD

echo "Retrieving Pristionchus genome assembly and annotation files"

#source 0README

echo "Performing alignment on the RNA-seq files"

readconfigfile $configfile

cd $genomedir

  if [ ! -e ${genomedir}/SAindex ] ; then
    echo "STAR --runMode genomeGenerate --runThreadN $numproc  ${STARgenomeGenerateOptions}  --genomeDir $genomedir --genomeFastaFiles $genomedir/*.fa"
    STAR --runMode genomeGenerate --runThreadN $nThreads  ${STARgenomeGenerateOptions}  --genomeDir $genomedir --genomeFastaFiles $genomedir/${genomeFasta}
  else
    echo "Using existing STAR suffix array for genome file $genomedir/*.fa"
  fi

  cd ${outDir}
  for file1 in ${WD}/${fqDir}/*.R1.fastq; do
    file2=$(basename $file1 .R1.fastq).R2.fastq
    echo "STAR --runMode alignReads --runThreadN $numproc  ${STARalignReadsOptions}  --outSAMtype BAM SortedByCoordinate --outSAMorder Paired  --outFileNamePrefix $(basename $file1 _001.fastq.gz).STAR.  --genomeDir $genomedir  --readFilesIn ${file1} ${file2}"
    STAR --runMode alignReads --runThreadN $nThreads  ${STARalignReadsOptions}  --outSAMtype BAM SortedByCoordinate --outSAMorder Paired  --outFileNamePrefix $(basename $file1 .fastq).STAR. --genomeDir $genomedir  --readFilesIn ${file1} ${file2}
  done

  cd ${WD}
  if [ ! -d "alignments" ] ; then
    mkdir alignments
  fi
  cd alignments
  for file1 in ${WD}/$fqDir/*.STAR.*.bam ; do
    echo $file1
    file2=$(basename $file1)
    echo $file2
    ln -s $file1 ./${file2/.STAR.Aligned.sortedByCoord.out}
  done
  cd ..

  echo ""
  echo " Done with step 3 (read mapping)."
  echo ""
  echo "================================================================================"
fi


exit


