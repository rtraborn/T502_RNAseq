
#creates a series of plots relating to differential expression analysis
# uses the bioconductor package 'limma'

## loading required packages (requires prior installation of each of the following)
require("limma") 
require("edgeR")
require("Rsubread")
require("Biobase")
require("gplots")
#require("DESeq2")

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

## end of read counts section ##

#load("prist_DE.RData") #starting from an R binary containing the featureCounts list created using the commands above. To run the above commands simply uncomment them (remove the leading '#' from each individual command), and commend out this line.

dge <- DGEList(counts = prist_fc$counts,
               group = c(rep("seud1",4),rep("nhr40",4)),
               genes = prist_fc$annotation$GeneID)
               
dge <- calcNormFactors(dge)

design <- model.matrix(~dge$samples$group)
               
colnames(design) <- c("seud1", "nhr40")

design #what does this object look like

dge <- estimateGLMCommonDisp(dge, design) #estimate the dispersion

dge <- estimateGLMTagwiseDisp(dge, design) 

plotBCV(dge)

#plotMDS(dge)

v <- voom(dge, design, plot=TRUE)

fit <- glmFit(dge, design) #fit the results to a linear model

lrt <- glmLRT(fit, coef = 2) #performs a likihood ratio test

save(dge, file="pristDGE.RData")

prist_top_tags <- topTags(lrt, n=1000, adjust.method="BH", sort.by="PValue", p.value=0.01)

save(prist_top_tags, file= "prist_top_tags.RData")

head(prist_top_tags[[1]]) #shows the top results on the screen

write.csv(prist_top_tags[[1]], file="prist_top_tags.csv", col.names=TRUE, row.names=FALSE) #writes a csv file to your working directory

dev.off()               
               
##############

#Creating an ExpressionSet object to perform the heatmap operations on               
#Mn_eset <-new("ExpressionSet", exprs=as.matrix(Dp_edger))
               
#de_data <- Dp_dge$pseudo.counts

#differential analysis results
#de_data <- cbind(de_data, prist_top_tags)

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
