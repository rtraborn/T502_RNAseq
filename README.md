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

Once it has passed through the queue the job should take 20-30 minutes to run on Carbonate. In the meantime you can proceed to the next part of the analysis.