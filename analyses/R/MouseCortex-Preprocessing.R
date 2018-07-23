library(SingleCellExperiment)
library(scater)

#read in data
sce <- readRDS("~/Downloads/2015-Zeisel-MouseCortex.RDS")
sce

example_sce <- calculateQCMetrics(sce)
colnames(colData(example_sce))
dim(example_sce)

#filter by features
keep_feature <- nexprs(example_sce, byrow=TRUE) >= 4
example_sce <- example_sce[keep_feature,]

#filter by cells
#keep.total <- isOutlier(example_sce$total_counts, nmads=.5, 
#                       type="lower", log=TRUE)
keep.total <- example_sce$total_counts > 1e3
keep.n <- example_sce$total_features_by_counts > 500
filtered <- example_sce[,keep.total & keep.n]
dim(filtered)

#normalize
sizeFactors(filtered) <- librarySizeFactors(filtered)
summary(sizeFactors(filtered))
norm <- normalize(filtered)
dim(norm)



#get counts and colData from preprocessed and normalized data
normlizedmousedata <- assay(norm, "logcounts") # count table with gene along the rows and cells along the columns
mouseMeta <- colData(norm)         # phenotypic table containing meta information about dataset

center <- colMeans(normlizedmousedata)
countTable <- normlizedmousedata - rep(center, rep.int(nrow(normlizedmousedata), ncol(normlizedmousedata)))

#save to .csv file
write.csv(countTable, file = "~/Desktop/2018summerintern/python/mouseCounts.csv")