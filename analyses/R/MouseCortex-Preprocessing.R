library(SingleCellExperiment)
library(scater)

#read in data
#sce <- readRDS("~/Downloads/2015-Zeisel-MouseCortex.RDS")
sce

example_sce <- calculateQCMetrics(sce)
colnames(colData(example_sce))
dim(example_sce)

#filter by features
keep_feature <- nexprs(example_sce, byrow=TRUE) >= 4
example_sce <- example_sce[keep_feature,]

#filter by cells
keep.total <- example_sce$total_counts > 1e4
keep.n <- example_sce$total_features_by_counts > 500
filtered <- example_sce[,keep.total & keep.n]
dim(filtered)

#normalize
sizeFactors(filtered) <- librarySizeFactors(filtered)
summary(sizeFactors(filtered))
norm <- normalize(filtered)
dim(norm)

#get counts and colData from preprocessed and normalized data
countTable <- assay(norm, "counts") # count table with gene along the rows and cells along the columns
mouseMeta <- colData(norm)         # phenotypic table containing meta information about dataset

#save to .csv file
write.csv(countTable, file = "mouseCounts.csv")
write.csv(mouseMeta, file = "mouseMeta.csv")
