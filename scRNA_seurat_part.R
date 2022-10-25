###Cell quality control
###Removal of ambient RNA using SoupX
#https://cellgeni.github.io/notebooks/html/new-10kPBMC-SoupX.html

library(Seurat)
library(ggplot2)
library(knitr)
library(dplyr)

library(SoupX)
library(DropletUtils)
library(DoubletFinder)

setwd("mypath")
getwd()

#read in h5 files (filtered and raw)
filt.matrix <- Read10X_h5("filtered_feature_bc_matrix.h5",use.names = T)
raw.matrix  <- Read10X_h5("raw_feature_bc_matrix.h5",use.names = T)
str(filt.matrix)
str(raw.matrix)

#create seurat obj
srat  <- CreateSeuratObject(counts = filt.matrix)
srat

#make soup obj
soup.channel  <- SoupChannel(raw.matrix, filt.matrix)
soup.channel

rm(filt.matrix)
rm(raw.matrix)
#make clusters
srat    <- SCTransform(srat, verbose = F) %>% RunPCA() %>% RunUMAP(dims = 1:40) %>% FindNeighbors(dims = 1:40)
srat    <- FindClusters(srat, verbose = T)

#add cluster info to the channel using setClusters
#setDR is useful for visualizations
meta    <- srat@meta.data
umap    <- srat@reductions$umap@cell.embeddings
soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)
head(meta)
head(soup.channel)

rm(srat)
rm(meta)
#calculating ambient RNA profile
soup.channel  <- autoEstCont(soup.channel)

#Genes with highest expression in background. 
#These are often enriched for ribosomal proteins.
#for check proteins
head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 20)

#We will use roundToInt option to make sure we output integer matrix.
adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)

#Finally, let’s write the directory with corrected read counts.
DropletUtils:::write10xCounts("soupX_E8RNA_filt", adj.matrix)

rm(soup.channel)

###create seurat obj###
adj.matrix <- Read10X("./soupX_E8RNA_filt")
srt.obj <- CreateSeuratObject(counts = adj.matrix, assay = "RNA", project = "E8RNA")
srt.obj
srt.obj@meta.data$Barcodes <- rownames(srt.obj@meta.data)

rm(adj.matrix)

#grep all the MT genes
#Calculate percentage of mitochondrial on ribosomal counts
# Percent of mitochondrial counts
grep(pattern ='(mypattern)',rownames(srt.obj@assays$RNA@counts),value = TRUE)

srt.obj[["percent.mt"]] <- PercentageFeatureSet(srt.obj, pattern = '(mypattern)')
str(srt.obj@meta.data)

# Percent of mitochondrial ribosomal
grep(pattern ='(mypattern)',rownames(srt.obj@assays$RNA@counts),value = TRUE)

srt.obj[["percent.ribo"]] <- PercentageFeatureSet(srt.obj, pattern = '(mypattern)')
str(srt.obj@meta.data)

#Filtering cells based on number of genes and transcripts detected
#Remove cells with too few gene detected or with too many UMI counts
#Set low and hight thresholds on the number of detected genes

RNA.max <- round(mean(srt.obj$nFeature_RNA) + 2 * sd(srt.obj$nFeature_RNA), digits = -2)
RNA.min <- round(mean(srt.obj$nFeature_RNA) - 2 * sd(srt.obj$nFeature_RNA), digits = -2)

# Set minimum parameters to 0 if negative value
if (RNA.min < 0){
  RNA.min <- 0
} else {
  RNA.min <- RNA.min
}

# Set hight threshold on the number of transcripts
Cell.QC.Stat <- srt.obj@meta.data
max.nCount_RNA.thr <- median(Cell.QC.Stat$nCount_RNA) + 3*mad(Cell.QC.Stat$nCount_RNA)

# Filter cells base on both metrics
srt.obj_subset <- subset(srt.obj, subset = nFeature_RNA < RNA.max & 
                           nFeature_RNA > RNA.min &
                           nCount_RNA < max.nCount_RNA.thr &  
                           percent.mt < 10)

srt.obj_subset
head(srt.obj_subset@meta.data)

###Use Scrublet to detect obvious doublets
#https://matthieuxmoreau.github.io/EarlyPallialNeurogenesis/html-Reports/Quality_Control.html#Use_Scrublet_to_detect_obvious_doublets
install.packages("ggExtra")

library(dplyr)
library(RColorBrewer)
library(ggExtra)
library(cowplot)
library(reticulate)

#Export filtered matrix to input matrix format for scrublet
library(DropletUtils)

write10xCounts(x = srt.obj_subset@assays$RNA@counts, path = '/home/miller/Desktop/ZF_scAnalysis/CellRanger_opt/scRNA/soupX_seuratQC_E8RNA_filt')

# Check the current Python version
# !!!scrublet has to be installed and used in the same Python version!!!
reticulate::py_config()
# mypath: /home/miller/.local/share/r-miniconda/envs/r-reticulate/bin

# Loads Python Shell
repl_python()

#(Python) Export raw count matrix as input to Scrublet
import scrublet as scr
import scipy.io
import numpy as np
import os

#(Python) Load raw counts matrix and gene list
input_dir = '/home/miller/Desktop/ZF_scAnalysis/CellRanger_opt/scRNA'
counts_matrix = scipy.io.mmread(input_dir + '/soupX_seuratQC_E8RNA_filt/matrix.mtx').T.tocsc()

#(Python) Initialize Scrublet object
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.1)

#(Python) Run the default pipeline
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

# Import scrublet's doublet score
srt.obj_subset@meta.data$Doubletscore <- py$doublet_scores

# Plot doublet score
ggplot(srt.obj_subset@meta.data, aes(x = Doubletscore, stat(ndensity))) +
  geom_histogram(bins = 200, colour ="lightgrey")+
  geom_vline(xintercept = 0.15, colour = "red", linetype = 2) # Manually set threshold

# Manually set threshold at doublet score to 0.15
srt.obj_subset@meta.data$Predicted_doublets <- ifelse(py$doublet_scores > 0.15, "Doublet","Singlet" )
table(srt.obj_subset@meta.data$Predicted_doublets)
head(srt.obj_subset@meta.data)

#remove doublet
srt.obj_subset_scrublet <- subset(srt.obj_subset, subset = Predicted_doublets=='Singlet')
