#!/bin/R

#creates a series of plots relating to differential expression analysis
# uses the bioconductor package 'limma'



require("limma") ## loading required packages (requires prior installation of each of the following)
require("edgeR")
require("Rsubread")
require("Biobase")
require("gplots")
require("DESeq2")

bamDir=/N/dc2/scratch/rtraborn/T502_fastq/PP_RNAseq
WD=/N/dc2/scratch/rtraborn/T502_RNAseq

#library #could be used in place of require()

#obtaining directory paths for both groups
seudDir <- "/scratch/rtraborn/archive/neuron_RNAseq/fastq/bam_alns/neuron_7mo" 
NHR40Dir <- c("/scratch/rtraborn/archive/neuron_RNAseq/fastq/bam_alns/neuron_13mo")

seud_files <- list.files(seudDir,pattern="\\.bam$", full.names=TRUE) #obtaining list of file names for both groups
NHR40_files <- list.files(NHR40Dir,pattern="\\.bam$", full.names=TRUE)

prist_files <- c(seud_files, NHR40_files)

#creating a count table
neuron_fc <- featureCounts(neuron_files, annot.inbuilt="mm10", useMetaFeatures=TRUE, strandSpecific=1, nthreads=6)

## end of read counts section ##

load("neuron_fc_obj.RData") #starting from an R binary containing the featureCounts list created using the commands above. To run the above commands simply uncomment them (remove the leading '#' from each individual command), and commend out this line.

dge <- DGEList(counts = neuron_fc$counts,
               group = c(rep("neuron_7",6),rep("neuron_13",6)),
               genes = neuron_fc$annotation$GeneID)
               
dge <- calcNormFactors(dge)

design <- model.matrix(~dge$samples$group)
               
colnames(design) <- c("neuron_7", "neuron_13")

design #what does this object look like

dge <- estimateGLMCommonDisp(dge, design) #estimate the dispersion

dge <- estimateGLMTagwiseDisp(dge, design) 

plotBCV(dge)

v <- voom(dge, design, plot=TRUE)

fit <- glmFit(dge, design) #fit the results to a linear model

lrt <- glmLRT(fit, coef = 2) #performs a likihood ratio test

neuro_top_tags <- topTags(lrt)

save(neuro_top_tags, file= "neuro_top_tags.RData")

head(neuro_top_tags) #shows the top results on the screen

write.csv(neuro_top_tags, file="neuro_top_tags.csv", col.names=TRUE, row.names=FALSE) #writes a csv file to your working directory

dev.off()               
               
##############

#Creating an ExpressionSet object to perform the heatmap operations on               
Mn_eset <-new("ExpressionSet", exprs=as.matrix(Dp_edger))
               
#de_data <- Dp_dge$pseudo.counts

#differential analysis results
#de_data <- cbind(de_data, neuro_top_tags)

#calculating the false discovery rate (FDR)
#de_data$FDR <- p.adjust(de_data$P.Value, method = 'BH')

#dispersion of each tag cluster
#de_data$tw_dis <- dge$tagwise.dispersion

#coordinates of each tag cluster
#data_coord2 <- matrix(data=unlist(strsplit(rownames(de_data), split="_")),
#                      nrow= length(row.names(de_data)),
#                      byrow=T)

#data_coord2 <- as.data.frame(data_coord2, stringsAsFactors=F)

#de_index <- which(de_data1$FDR<0.01)

#length(de_index)

#de <- de_table1[de_index,]

#selected <- rownames(de)

#esetSel <- Mn_eset[selected, ]

#heatmap.2(exprs(esetSel), symm=FALSE,symkey=FALSE,scale="row", density.info="none",trace="none",
#          key=TRUE,margins=c(10,10))

#dev.off()
