library(SingleCellExperiment)
library(scater)
library(scran)

#read in data
sce <- readRDS("data/2015-Zeisel-MouseCortex.RDS")
sce

ercc_genes <- grep("ERCC", rownames(sce))
mt_genes <- grep("mt-", rownames(sce))
sce <- calculateQCMetrics(sce, 
                  feature_controls = c(ERCC = ercc_genes, mt = mt_genes))
dim(sce)
table(sce$level1class) # smallest group has 98 cells

#### filter by features

# Here we want to look for genes that have least one count in 98 cells (smallest group identified by Zeisel et al)
keep_feature <- nexprs(sce, byrow=TRUE) >= 98
# Next we remove lowly expressed genes using a mean-based filter
ave.counts <- rowMeans(counts(sce))
keep_feature <- keep_feature & rowMeans(counts(sce)) >= 0.2
sum(keep_feature) # we will keep 8912 genes

#### filter by cells
# libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
# feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)

sce$use <- (sce$total_features > 1000 & # sufficient features (genes)
            sce$total_counts > 5000 # sufficient molecules counted
)
table(sce$use) # keep 2992 cells
filtered <- sce[keep_feature, sce$use]
dim(filtered)

#### normalize
## size factor normalisation with scran 
qclust <- quickCluster(filtered) # this step is slow. It is also important because the cell types are very, very different. The assumptions fo the normalization require it. 
filtered <- computeSumFactors(filtered, clusters = qclust)
summary(sizeFactors(filtered))

# wider scatter plot than just along the line; reflects DE between cell types in Zeisel data
plot(sizeFactors(filtered), filtered$total_counts/1e3, log="xy",
     ylab="Library size (thousands)", xlab="Size factor")

# compute scaling factors (total counts of ERCC spikes) separately for ERCC spikes
filtered <- computeSpikeFactors(filtered, type="ERCC", general.use=FALSE)

# normalize data
filtered <- normalize(filtered)

#### remove ERCC spike in and MT genes for downstream analyses
filtered <- filtered[!rowData(filtered)$is_feature_control, ]

# test out tSNE in R
plotTSNE(filtered, colour_by = "level1class", rand_seed = 20180808)
# see quite a bit of clustering in the tSNE plot

#### Save data
#get counts and colData from preprocessed and normalized data
normlizedmousedata <- assay(filtered, "logcounts") # count table with gene along the rows and cells along the columns
mouseMeta <- colData(filtered)         # phenotypic table containing meta information about dataset

center <- rowMeans(normlizedmousedata)
countTable <- sweep(normlizedmousedata, 1, center, FUN="-")
#countTable <- normlizedmousedata - rep(center, rep.int(nrow(normlizedmousedata), ncol(normlizedmousedata)))

#save to .csv file
write.csv(countTable, file = "~/Desktop/2018summerintern/python/mouseData.csv")
write.csv(colData(filtered)["level1class"], file = "~/Desktop/2018summerintern/python/realMouse.csv")


# Useful tutorials for Zeisel data set
# https://github.com/davismcc/scater_tutorials_open_data/blob/master/zeisel_mouse_cortex.Rmd
# https://f1000research.com/articles/5-2122/v2