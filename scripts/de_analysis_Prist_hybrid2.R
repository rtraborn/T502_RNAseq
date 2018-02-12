
#creates a series of plots relating to differential expression analysis
# uses the bioconductor package 'limma'

## loading required packages (requires prior installation of each of the following)
require("limma") 
require("edgeR")
require("Rsubread")
require("Biobase")
require("gplots")
require("mixomics")
#require("DESeq2")


### set your path to the scripts directory
WD="/N/dc2/scratch/rtraborn/T502_RNAseq/scripts"
setwd(WD)

#obtaining directory paths for both groups
bamDir <- "../alignments"
pristAnnot <- "../annotation/Hybrid2_AUGUSTUS2014_gene.gtf"

#obtaining list of file names for both groups
seud_files <- list.files(bamDir, pattern="\\Seud1", full.names=TRUE)
NHR40_files <- list.files(bamDir, pattern="\\NHR40", full.names=TRUE)

prist_files <- c(seud_files, NHR40_files)

#creating a count table
prist_fc <- featureCounts(prist_files, annot.ext=pristAnnot, useMetaFeatures=TRUE, strandSpecific=1, isPairedEnd=FALSE, nthreads=16, isGTFAnnotationFile=TRUE, primaryOnly=TRUE)

save(prist_fc, file="prist_DE.RData") #saving our featureCounts data to an R binary

## end of read counts section ##

#load("prist_DE.RData") #starting from an R binary containing the featureCounts list created using the commands above. To run the above commands simply uncomment them (remove the leading '#' from each individual command), and commend out this line.

dge <- DGEList(counts = prist_fc$counts,
               group = c(rep("seud1",4),rep("nhr40",4)),
               genes = prist_fc$annotation$GeneID)

### Now we will apply TMM normalization
               
dge <- calcNormFactors(dge)

### Let's take a look at the normalization factors

dge$samples

# Filtering out genes that have a low number of counts (i.e. are lowly-expressed)

keep <- rowSums(cpm(dge)>1) >= 2
dge <- dge[keep, keep.lib.sizes=FALSE]

# Creating a design matrix to model our experiment

design <- model.matrix(~dge$samples$group)

colnames(design) <- c("seud1", "nhr40")

design #what does this object look like?

#estimate the dispersion

dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design) 

# evaluate the common dispersion

sqrt(dge$common.disp)

# Now we plot the tagwise dispersions against the log2-scaled counts-per million (CPM) values

plotBCV(dge)

# Now we performed the differential expression calculation

fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

summary(de <- decideTestsDGE(lrt, p=0.01, adjust="BH"))
de_tags <- rownames(decideTestsDGE(lrt, p=0.01, adjust="BH"))
de_tags <- rownames(dge)[as.logical(de)]

#mkaing a smear (i.e. a mean-difference) plot of our data

plotSmear(lrt, de.tags=de_tags)
abline(h=c(-2,2), col="blue")

# We can also make the (classic) volcano plot from our data
volcanoData <- cbind(lrt$table$logFC, -log10(lrt$table$PValue))
plot(volcanoData, pch=19)
abline(v=c(-2,2), col="red")

save(dge, file="pristDGE.RData") #saving the updated dge object to our working directory

prist_top_tags <- topTags(lrt, adjust.method="BH", sort.by="PValue", p.value=0.01)
head(prist_top_tags[[1]]) #shows the top results on the screen
write.csv(prist_top_tags[[1]], file="prist_top_tags.csv", row.names=FALSE) #writes a csv file to your working directory
save(prist_top_tags, file= "prist_top_tags.RData") #saves the prist_top_tags file as a p-value

##############
#Making a heatmap with the differentially-expressed genes
library(Biobase) #load this required package if you haven't already done so
library(gplots)

de_data <- dge$counts
colnames(de_data) <- c("Seud1-1","Seud1-2","Seud1-3", "Seud1-4", "NHR40-1","NHR40-2", "NHR40-3", "NHR40-4")
head(de_data)

top_tags <- topTags(lrt, n= 18146, sort.by="none")

#differential analysis results
de_data <- cbind(de_data, top_tags)

#calculating the false discovery rate (FDR)
de_data$FDR <- p.adjust(de_data$P.Value, method = 'BH')

diff.genes = rownames(de_data[de_data$FDR<0.01, ])
head(diff.genes)
length(diff.genes)

dge.subset = dge[diff.genes, ]

heatmap.2(dge.subset$counts, symm=FALSE,symkey=FALSE,scale="row", density.info="none",trace="none",
          key=TRUE,margins=c(10,10))

