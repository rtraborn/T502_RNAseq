# RNA-seq workflow for BIOT-T502

## Instructions

In today's class (Tuesday the 30th) we will run the first part of our RNA-seq analysis pipeline. To do this we will use the many of the skills that we learned over the previous three class periods. 

To begin, please click on this repository on the github: https://github.com/rtraborn/T502_RNAseq . Then, please click 'Fork' on the upper-right-hand side of the page. This will fork this github repo and it attach it to your own personal github account, allowing you to make changes and push them to your own account. 

Now, please navigate to your own personal github page. Under 'repositories', you should see _T502_RNA-seq_. Please click on it. On the main repository page you should see *<your user ID>*/T502_RNAseq in the upper-left-hand corner of the page, and underneath it, in smaller letters, you should see "forked from rtraborn/T502_RNA-seq" in smaller letters.

Once that checks out we can go head and clone the repository to your home directory on Carbonate. Please use your ssh2 client to login to Carbonate. Then, once you have logged in, please copy the link to your forked repository. It should look something like this: https://github.com/studentID/T502_RNA-seq (where _studentID_ is your github account username).

Once logged in, now clone the repository to your home directory in Carbonate, and enter the directory as follows:

```
git clone https://github.com/studentID/T502_RNA-seq
cd T502_RNA-seq 
```

Now, let's enter the `scripts/` directory: 

```
cd scripts
```

### Running Fastqc on our reads

If type `ls` in our bash interpreter we'll see that there are several scripts in this directory, but let's focus our attention on the one entitled `fastqc_batch.sh`. We can read the contents of the script using by typing `less fastqc_batch.sh`, which we can exit at any time by typing the letter 'q'. As you'll see, this script will make links to the RNAseq files hosted on Carbonate, and then will complete a [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) run on all of them.

It will be necessary to make a few changes to our repository before we can submit this job to the scheduler. After the PBS header lines, we'll need to specify a output directory.

```
####### Before running the script, please enter path to desired output directory, below ####  
fqDir=<provide path to desired output directory>  
```

Please type `nano fastqc_batch.sh` and use the text editor to enter your desired output directory on Carbonate. The script will not run until you've made this change.

Once you've done this we can submit our Fastqc job to the scheduler.

```
qsub fastqc_batch.sh
```

Once it has passed through the queue the job should take 20-30 minutes to run on Carbonate. Once this is complete please download the results (which will be in .zip format) to your local workstation. Once they are on your computer, uncompress them and double-click on the *.html file inside the folder to view each report. We will go through how to interpret these results together during class.

### Alignment

Now let's proceed to the alignment, which will be the most computatationally-expensive job we will run in our analysis. We will be taking both reads pairs from all eight (8) isolates and align them to the _Pristionchus pacificus_ genome. Of course, this will require that we also download the genome assembly. I have put together everything we'll need for this job, including a script (found in `./0README`) that automatically downloads the assembly and annotation. The aligner software we will use is called STAR (Dobin et al., 2013). The documentation for STAR can be found here: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf. I have chosen to run STAR because it is fast, accurate and highly customizable. Fortunately, this software is installed as a module in Carbonate, so we can invoke it easily in our batch script (`module load star` to be precise) without having to install it ourselves.

Before we begin, make sure that you've commited the changes that you've already made in your repository. It's important that you stage (add), commit and push all of your updated files (especially the fastqc script) before proceeding. Assuming you are in this directory, please type:

```
git status #see list of updated files
git add <path to file(s) to add> #replace everything after "add" with path(s) to the file(s) in your repo that you want to track (*not* fastq files!). 
git commit -m "Compose a descriptive commit message."
put push origin master
```

Great. Now that this is complete you'll need to pull from this (i.e. my) repository, as I've added a few important updates to the workflow (and to this document also!). 

``` 
## assuming you are still in the workflow directory:
git pull https://github.com/rtraborn/T502_RNAseq master
```

It may try to do a 'merge' on some of your recently-changed files. If it does this you'll need to select which of the two versions of a file to keep. In the case of the `fastqc_batch.sh` script, choose your changes. Everything else should go through ok, and it will sync my recent updates. 

First, let's automatically download the assembly and annotation files.

```
source 0README
## this will initiate download of the necessary files to this repository.
ls 
```

Next, please move into the `scripts/` directory and edit the file using your favorite text editor.

```
## assuming you use nano
nano align_STAR.sh
```

Please make the necessary changes to the following lines:
```
WD=/N/u/<yourUserId>/Carbonate/T502_RNAseq/ 
```

There are are few things to note about the script (`align_STAR.sh`) that you are about to run. First, the script uses a configuration file `STARalign.conf` that contains the parameters we'll use for our STAR alignment. Second, if you look at the header of the script you'll see that the wall time (i.e. the amount of time we're requesting) is higher than before. As you might expect, this is because of the length of time it takes to run the alignment. I did mention that STAR is fast, but the genome indexing portion of the run (the first STAR command in the script) takes a pretty decent amount of time.

Once we've checked over the script and added the correct path to the line beginning in `WD=`, then we can submit our job to the scheduler.

```
#assuming we're in scripts/
qsub align_STAR.sh
```

Log back into Carbonate periodically to check on the status of your job as follows:

```
qstat -u yourID
```

If a job terminates prematuresly, please look at the error file (e.g. `PP_RNAseq_STAR_align_T502.e109469`) in your scripts directory, which will give you clues as to why the job halted.

```
### your error file name will be slighly different- this is just an example.
less PP_RNAseq_STAR_align_T502.e109469 
```

Once your job is complete, you should have 8 BAM files in the `alignments/` directory. You can look at them using a module called samtools.

```
cd alignments
module load samtools
samtools view NHR40-2.R1.bam | less
## that should open a window in your terminal that allows you to explore the alignments using your cursor.
## exit by typing the letter 'q' at any time.
```

### Differential gene expression analysis

Now we will analyze our data to determine the set of genes that are upregulated. To do this, we will need to integrate our alignment data with the set of gene features from Pristionchus.

The repo has been updated since the previous class. To bring your repo up to date, please stage, commit and push your changes (not log files, alignment files or fastqs, though!) to your fork. See above for instructions.

Next, you'll need to pull from my repo as follows:
```
git pull https://github.com/rtraborn/T502_RNAseq master
### type :wq to exit the (vim) merge screen

```


#### Installation instructions

Before we can continue, it is necessary to install a set of R packages, most of them from the Bioconductor package repository. Do do this, please following the instructions in `./INSTALL.md` by typing `less INSTALL.md` from the top level of this directory. It will be necessary to do this from the R console in both RStudio and Carbonate.

#### Submitting the job on Carbonate
If you want to perform the entireity of this job on Carbonate, we will use the batch script, `de_launch.sh` that is included in `scripts/`.
Similarly, you will also need to edit the 14th line of the R script `de_analysis_Prist_hybrid2.R`.
Please edit the file using a text editor (as described previously) to provide the necessary path on your file system.

Once this is complete, please type the following:
```
cd scripts
qsub de_launch.sh
```

During class we will go through the analysis, and (explore the results together) using RStudio.

#### Interactive DGE analysis using an R console in RStudio

The following commands will help get you started:

```
#loading the library files:
library("limma") 
library("edgeR")
library("Rsubread")
library("Biobase")
library("gplots")
```

Now we will load the R binary file. This contains the count data calculated from the bam alignments. Since it takes a bit of time, we will use the caculated counts matrix directly.

```
load("prist_DE.RData")
```