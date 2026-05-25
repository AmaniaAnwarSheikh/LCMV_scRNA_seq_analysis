# =========================================================
# Workflow B Style scRNA-seq Pipeline
# Arm_d8 vs Arm_d21 vs Clone13_d21
# Seurat v5 + SCTransform Workflow
# =========================================================


# =========================================================
# Install CRAN packages
# =========================================================

install.packages(c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "patchwork",
  "tidyverse",
  "BiocManager"
))
install.packages("RcppEigen")
install.packages("RSpectra")
install.packages("stringi")
install.packages("hdf5r")


install.packages('BiocManager')
BiocManager::install('glmGamPoi')
install.packages("sctransform")
# =========================================================
# Install Bioconductor packages
# =========================================================

BiocManager::install(c(
  "scran",
  "org.Mm.eg.db",
  "SingleCellExperiment"
))

a
# =========================================================
# Load packages
# =========================================================

library(Seurat)
library(sctransform)
library(uwot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(scran)
library(org.Mm.eg.db)
library(SingleCellExperiment)


sessionInfo()

# =========================================================
# Create output folders
# =========================================================

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("objects", showWarnings = FALSE)

# =========================================================
# Define data directory
# =========================================================

data_dir <- "/Users/sheikh.a/Desktop/scRNA data"

# =========================================================
# Load Batch 1 files
# =========================================================

Arm_d8_B1.data <- Read10X_h5(
  file.path(
    data_dir,
    "GSM7919089_SC_Batch_1_Arm_d8_lane_A_sample_feature_bc_matrix.h5"
  )
)

Arm_d21_B1.data <- Read10X_h5(
  file.path(
    data_dir,
    "GSM7061752_SC_Batch_1_Arm_d21_lane_Asample_feature_bc_matrix.h5"
  )
)

Clone13_d21_B1.data <- Read10X_h5(
  file.path(
    data_dir,
    "GSM7061754_SC_Batch_1_Clone13_d21_lane_Asample_feature_bc_matrix.h5"
  )
)

# =========================================================
# Load Batch 2 files
# =========================================================

Arm_d8_B2.data <- Read10X_h5(
  file.path(
    data_dir,
    "GSM7919090_SC_Batch_2_Arm_d8_sample_feature_bc_matrix.h5"
  )
)

Arm_d21_B2.data <- Read10X_h5(
  file.path(
    data_dir,
    "GSM7061756_SC_Batch_2_Arm_d21_sample_feature_bc_matrix.h5"
  )
)

Clone13_d21_B2.data <- Read10X_h5(
  file.path(
    data_dir,
    "GSM7061757_SC_Batch_2_Clone13_d21_sample_feature_bc_matrix.h5"
  )
)

# =========================================================
# Create Seurat objects
# =========================================================

Arm_d8_B1 <- CreateSeuratObject(
  counts = Arm_d8_B1.data$`Gene Expression`,
  project = "Arm_d8"
)

Arm_d21_B1 <- CreateSeuratObject(
  counts = Arm_d21_B1.data$`Gene Expression`,
  project = "Arm_d21"
)

Clone13_d21_B1 <- CreateSeuratObject(
  counts = Clone13_d21_B1.data$`Gene Expression`,
  project = "Clone13_d21"
)

Arm_d8_B2 <- CreateSeuratObject(
  counts = Arm_d8_B2.data$`Gene Expression`,
  project = "Arm_d8"
)

Arm_d21_B2 <- CreateSeuratObject(
  counts = Arm_d21_B2.data$`Gene Expression`,
  project = "Arm_d21"
)

Clone13_d21_B2 <- CreateSeuratObject(
  counts = Clone13_d21_B2.data$`Gene Expression`,
  project = "Clone13_d21"
)

# =========================================================
# Add metadata
# =========================================================

# Arm_d8
Arm_d8_B1$condition <- "Arm"
Arm_d8_B1$day <- 8
Arm_d8_B1$batch <- "Batch1"

Arm_d8_B2$condition <- "Arm"
Arm_d8_B2$day <- 8
Arm_d8_B2$batch <- "Batch2"

# Arm_d21
Arm_d21_B1$condition <- "Arm"
Arm_d21_B1$day <- 21
Arm_d21_B1$batch <- "Batch1"

Arm_d21_B2$condition <- "Arm"
Arm_d21_B2$day <- 21
Arm_d21_B2$batch <- "Batch2"

# Clone13_d21
Clone13_d21_B1$condition <- "Clone13"
Clone13_d21_B1$day <- 21
Clone13_d21_B1$batch <- "Batch1"

Clone13_d21_B2$condition <- "Clone13"
Clone13_d21_B2$day <- 21
Clone13_d21_B2$batch <- "Batch2"

# =========================================================
# Merge all datasets
# =========================================================

combined <- merge(
  Arm_d8_B1,
  y = list(
    Arm_d8_B2,
    Arm_d21_B1,
    Arm_d21_B2,
    Clone13_d21_B1,
    Clone13_d21_B2
  ),
  add.cell.ids = c(
    "Arm_d8_B1",
    "Arm_d8_B2",
    "Arm_d21_B1",
    "Arm_d21_B2",
    "Clone13_B1",
    "Clone13_B2"
  )
)

# =========================================================
# Calculate mitochondrial percentage
# =========================================================

combined[["percent.mt"]] <- PercentageFeatureSet(
  combined,
  pattern = "^mt-|^MT-"
)

# =========================================================
# QC violin plots
# =========================================================

VlnPlot(
  combined,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  split.by = "batch",
  ncol = 3,
  pt.size = 0.1
)

ggsave(
  "figures/QC_violin_plots.png",
  width = 14,
  height = 5,
  dpi = 300
)

# =========================================================
# Filter cells
# =========================================================

combined <- subset(
  combined,
  subset =
    nFeature_RNA > 800 &
    nFeature_RNA < 6000 &
    percent.mt < 7
)

# =========================================================
# SCTransform normalization
# =========================================================

combined <- SCTransform(
  combined,
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

# =========================================================
# Remove TCR genes from variable features
# =========================================================

varfeatures <- combined@assays$SCT@var.features

varfeatures <- varfeatures[
  !grepl("^Tr[ab]", varfeatures)
]

combined@assays$SCT@var.features <- varfeatures

# =========================================================
# PCA
# =========================================================

combined <- RunPCA(
  combined,
  features = varfeatures,
  npcs = 50,
  verbose = FALSE
)

# =========================================================
# Elbow plot
# =========================================================

ElbowPlot(combined, ndims = 50)

ggsave(
  "figures/ElbowPlot.png",
  width = 6,
  height = 4,
  dpi = 300
)

# =========================================================
# UMAP and clustering
# =========================================================

combined <- RunUMAP(
  combined,
  dims = 1:25,
  verbose = FALSE
)

combined <- FindNeighbors(
  combined,
  dims = 1:25,
  verbose = FALSE
)

combined <- FindClusters(
  combined,
  resolution = 0.5,
  verbose = FALSE
)

# =========================================================
# UMAP visualization
# =========================================================

p1 <- DimPlot(
  combined,
  group.by = "orig.ident"
) + ggtitle("By Sample")

p2 <- DimPlot(
  combined,
  group.by = "batch"
) + ggtitle("By Batch")

p3 <- DimPlot(
  combined,
  label = TRUE
) + ggtitle("By Cluster")

(p1 | p2 | p3)

ggsave(
  "figures/Combined_UMAPs.png",
  width = 18,
  height = 5,
  dpi = 300
)

# =========================================================
# Cell cycle analysis
# =========================================================

mm.pairs <- readRDS(
  system.file(
    "exdata",
    "mouse_cycle_markers.rds",
    package = "scran"
  )
)

mm.pairs.symbol <- lapply(
  mm.pairs,
  function(x) {
    
    data.frame(
      first = mapIds(
        org.Mm.eg.db,
        x$first,
        "ENSEMBL",
        column = "SYMBOL"
      ),
      second = mapIds(
        org.Mm.eg.db,
        x$second,
        "ENSEMBL",
        column = "SYMBOL"
      )
    )
  }
)

# =========================================================
# Convert Seurat object
# =========================================================

combined_sce <- as.SingleCellExperiment(
  combined,
  assay = "SCT"
)

# =========================================================
# Run cyclone
# =========================================================

library(scran)
library(SingleCellExperiment)

# Extract matrix
mat <- assay(combined_sce, "logcounts")

# Filter genes
keep <- rowSums(mat > 0) >= 20
mat <- mat[keep, ]

# Keep only cyclone genes
cyclone.genes <- unique(unlist(mm.pairs.symbol))
mat_small <- mat[rownames(mat) %in% cyclone.genes, ]

dim(mat_small)

# Run cyclone WITHOUT parallelization
assigned <- scran::cyclone(
  mat_small,
  pairs = mm.pairs.symbol
)

# Store phases
combined_sce$Phase <- assigned$phases

# Check results
table(combined_sce$Phase)

# Save
saveRDS(assigned, "cyclone_assignment.rds")
saveRDS(combined_sce, "combined_sce_with_cycle.rds")

# =========================================================
# Add phase metadata
# =========================================================

combined$Phase <- assigned$phases
# =========================================================
# Visualize cell cycle phases
# =========================================================

DimPlot(
  combined,
  group.by = "Phase"
)

ggsave(
  "figures/CellCycle_Phases.png",
  width = 6,
  height = 5,
  dpi = 300
)


# =========================================================
# Marker analysis
# =========================================================

# Set cluster identities
Idents(combined) <- "seurat_clusters"

# IMPORTANT:
# Required for Seurat v5 SCT differential expression
combined <- PrepSCTFindMarkers(combined)

# Find markers
markers <- FindAllMarkers(
  combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Check marker table
dim(markers)
head(markers)

# Save all markers
write.csv(
  markers,
  "results/cluster_markers.csv",
  row.names = FALSE
)

# =========================================================
# Top markers per cluster
# =========================================================

library(dplyr)

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  )

# Save top markers
write.csv(
  top10,
  "results/top10_markers.csv",
  row.names = FALSE
)

# View top markers
head(top10)


# =========================================================
# T-cell marker visualization
# =========================================================

genes <- c(
  "Klrg1","Zeb2","Cx3cr1","Gzmb","Prf1",
  "Il7r","Ccr7","Sell","Tcf7","Klf2",
  "Cd69","Itgae","Cxcr6",
  "Bcl6","Cxcr5","Pdcd1","Icos","Il21",
  "Havcr2","Lag3","Tigit","Ctla4","Tox"
)

genes <- genes[
  genes %in% rownames(combined)
]

DotPlot(
  combined,
  features = genes
) + RotatedAxis()

ggsave(
  "figures/Tcell_marker_dotplot.png",
  width = 12,
  height = 6,
  dpi = 300
)

# =========================================================
# OPTIONAL:
# Cell cycle regression
# =========================================================

# =========================================================
# OPTIONAL:
# Cell cycle regression
# =========================================================

library(Seurat)
library(sctransform)
library(future)

# IMPORTANT:
# Prevent future/memory errors
future::plan("sequential")

options(
  future.globals.maxSize = 20 * 1024^3
)

# Add cell cycle scores
combined$S_score <- assigned$scores$S
combined$G2M_score <- assigned$scores$G2M

# Regress cell cycle effects
combined_cc <- SCTransform(
  combined,
  vars.to.regress = c(
    "S_score",
    "G2M_score"
  ),
  verbose = FALSE
)

# =========================================================
# PCA after regression
# =========================================================

varfeatures_cc <- combined_cc@assays$SCT@var.features

varfeatures_cc <- varfeatures_cc[
  !grepl("^Tr[ab]", varfeatures_cc)
]

combined_cc@assays$SCT@var.features <- varfeatures_cc

combined_cc <- RunPCA(
  combined_cc,
  features = varfeatures_cc,
  npcs = 50,
  verbose = FALSE
)

combined_cc <- RunUMAP(
  combined_cc,
  dims = 1:25,
  verbose = FALSE
)

combined_cc <- FindNeighbors(
  combined_cc,
  dims = 1:25,
  verbose = FALSE
)

combined_cc <- FindClusters(
  combined_cc,
  resolution = 0.5,
  verbose = FALSE
)

# =========================================================
# Compare clustering after regression
# =========================================================

DimPlot(
  combined_cc,
  group.by = c(
    "orig.ident",
    "batch",
    "Phase",
    "seurat_clusters"
  ),
  label = TRUE
)

ggsave(
  "figures/CellCycle_Regressed_UMAP.png",
  width = 16,
  height = 12,
  dpi = 300
)

# =========================================================
# Manual annotation
# =========================================================

map <- c(
  "0" = "Effector T cells",
  "1" = "Exhausted T cells",
  "2" = "Memory T cells",
  "3" = "Tpex cells"
)

lev <- levels(Idents(combined))
new_lev <- lev

for (k in names(map)) {
  if (k %in% lev) {
    new_lev[lev == k] <- map[[k]]
  }
}

levels(Idents(combined)) <- new_lev

combined$cluster_label <- Idents(combined)

# =========================================================
# Final annotated UMAP
# =========================================================

final_umap <- DimPlot(
  combined,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) + ggtitle("Annotated T-cell States")

print(final_umap)

ggsave(
  "figures/Annotated_UMAP.png",
  width = 8,
  height = 6,
  dpi = 300
)

# =========================================================
# Save Seurat objects
# =========================================================

saveRDS(
  combined,
  file = "objects/combined_SCTransform.rds"
)

saveRDS(
  combined_cc,
  file = "objects/combined_SCTransform_CCregressed.rds"
)



dir.create("objects", showWarnings = FALSE)

saveRDS(
  combined,
  file = "objects/combined_SCTransform.rds"
)




# =========================================================
# Finished
# =========================================================

message("Workflow complete.")


library(Seurat)
library(ShinyCell)

# Reload Seurat object
combined <- readRDS(
  "objects/combined_SCTransform.rds"
)

# Use SCT assay
DefaultAssay(combined) <- "SCT"

# Create lightweight object
combined_shiny <- DietSeurat(
  combined,
  assays = "SCT",
  dimreducs = c("pca", "umap"),
  graphs = NULL
)

# Set identities
Idents(combined_shiny) <- "cluster_label"

# Create config
scConf <- createConfig(combined_shiny)

# Generate app
makeShinyApp(
  combined_shiny,
  scConf,
  gene.mapping = FALSE,
  gex.assay = "SCT",
  shiny.title = "LCMV CD8 T-cell Atlas",
  shiny.dir = "LCMV_ShinyCell"
)

install.packages("systemfonts")
install.packages("textshaping")
install.packages(c("ragg",
  "svglite",
  "Cairo"
))

shiny::runApp("LCMV_ShinyCell")
