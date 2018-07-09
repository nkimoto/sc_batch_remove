# Usage  : Rscript Seurat_Combat.R data_batch1 data_batch2 ...
# Author : kimoton

library(Seurat)
library(dplyr)
library(Matrix)
library(sva)


# Input Datas
args <- commandArgs(trailingOnly = TRUE)
filename <- c()
filename <- append(filename, args[1])
data1 <- Read10X(data.dir = args[1])
sobj1 <- CreateSeuratObject(raw.data=data1,
                            project=filename[1])
data.combind <- sobj1
for (i in 2:length(args)){
    assign(paste0("data", i), Read10X(data.dir = args[i]))
    filename <- append(filename, args[i])
    data.combind <- AddSamples(object = data.combind,
                               new.data = get(paste0("data", i)),
                               add.cell.id = filename[i])
}
bp <- data.combind

cat("\n========== batch distribution before filtering ===========")
table(bp@meta.data$orig.ident)
mito.genes <- grep(pattern = "^MT-",
                   x = rownames(x = bp@data),
                   value = TRUE)
percent.mito <- Matrix::colSums(bp@raw.data[mito.genes, ]) / Matrix::colSums(bp@data)

# AddMetaData adds columns to object@data.info
bp <- AddMetaData(object = bp, metadata = percent.mito, col.name = "percent.mito")

# Filter by added meta data
bp <- FilterCells(object = bp,
                  subset.names = c("percent.mito"),
                  low.thresholds = c(-Inf),
                  high.thresholds = c(0.05))
bp <- FilterCells(object = bp,
                  subset.names = c("nGene"),
                  low.thresholds = c(500),
                  high.thresholds = c(2000))

# Normalize
bp <- NormalizeData(object = bp,
                    normalization.method = "LogNormalize",
                    scale.factor = 1e4)
## Filter by variance
bp <- FindVariableGenes(object = bp,
                        mean.function = ExpMean,
                        dispersion.function = LogVMR,
                        x.low.cutoff = 0.2,
                        x.high.cutoff = 4,
                        y.cutoff = 1.0)

length(x = bp@var.genes)

bp <- ScaleData(object = bp)

bp <- RunPCA(object = bp,
             pc.genes = bp@var.genes,
             do.print = FALSE,
             pcs.print = 1:2,
             pcs.compute = 40,
             maxit = 500,
             weight.by.var = FALSE)
pdf("pca.pdf")
PCAPlot(object = bp, dim.1 = 1, dim.2 = 2)
dev.off()

# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
bp <- ProjectPCA(object = bp, do.print = FALSE)

bp <- seuset <- JackStraw(object = bp,
                          num.replicate = 100)
pdf("jackstraw.pdf")
JackStrawPlot(object = bp, PCs = 1:20)
dev.off()

## Execute clustering
bp <- FindClusters(object = bp,
                   reduction.type = "pca",
                   dims.use = 1:10,
                   resolution = 0.6,
                   print.output = 0,
                   save.SNN = TRUE)

## Print cluster params
PrintFindClustersParams(object = bp)

## Run tSNE
bp <- RunTSNE(object = bp, dims.use = 1:10, do.fast = TRUE)
pdf("tsne.pdf")
TSNEPlot(object = bp)
dev.off()

### Run tSNE marked by nUMI
pdf("tsne_numi.pdf")
FeaturePlot(bp, features.plot=c('nUMI'))
dev.off()


batchname = bp@meta.data$orig.ident
batchid = rep(1,length(batchname))
batchid[batchname=="pbmc_4k_filtered"] = 2
names(batchid) = rownames(bp@meta.data)
bp <- AddMetaData(object = bp,
                  metadata = batchid,
                  col.name = "batchid")

cat("\n========== batch distribution after filtering ===========")
table(bp@meta.data$batchid)
cat("\n")


pdf("tsne_before_correct.pdf")
FeaturePlot(bp, features.plot=c('batchid'))
dev.off()
save(bp, file = "bp_pre_batch_correct.Robj")


load("bp_pre_batch_correct.Robj")

# after filtering cells, some genes have zero counts
bp@data = bp@data[Matrix::rowSums(bp@data)>0,] 
bp <- ScaleData(object = bp,
                vars.to.regress = c("batchid"))


## SEURAT
bp <- RunPCA(object = bp,
             pc.genes = bp@var.genes,
             do.print = FALSE,
             pcs.compute = 40,
             weight.by.var = FALSE)
bp <- RunTSNE(object = bp,
              dims.use = 1:10,
              do.fast = TRUE)

pdf("tsne_seurat.pdf")
FeaturePlot(bp, features.plot=c('batchid'))
dev.off()

## COMBAT
## restore it before running combat on it

load("bp_pre_batch_correct.Robj")
m = as.data.frame(as.matrix(bp@data))
m = m[rowSums(m)>0,]
com = ComBat(as.matrix(m),
             batchid,
             prior.plots=FALSE,
             par.prior=TRUE)

## Standardizing
bp@data = Matrix(as.matrix(com))
bp = ScaleData(bp)

bp <- RunPCA(object = bp,
             pc.genes = bp@var.genes,
             do.print = FALSE,
             pcs.print = 1:2,
             pcs.compute = 40,
             weight.by.var = FALSE)
bp <- RunTSNE(object = bp,
              dims.use = 1:10,
              do.fast = TRUE)

pdf("tsne_combat.pdf")
FeaturePlot(bp, features.plot=c('batchid'))
dev.off()
