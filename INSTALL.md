## Installing R packages required for the differential gene expression analysis.

To get an R console in Carbonate, please type `module load r' followed by `R`.
From an R console, please type the following:

### Installing the Bioconductor packages:
```
source("https://bioconductor.org/biocLite.R")
biocLite("limma") #installing limma
biocLite("edgeR") #installing edgeR
biocLite("Rsubread") #installing Rsubread
biocLite("Biobase") #installing Biobase
```

### now, the CRAN package we need:
```
install.packages("gplots")
```
