
# =========================================================
# scRNA-seq Analysis Pipeline
# LCMV Arm vs Clone13
# Compatible with Seurat v5
# =========================================================

# =========================================================
# Load packages
# =========================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

sessionInfo()

# =========================================================
# Define raw data location
# =========================================================

data_dir <- "/Users/sheikh.a/Desktop/scRNA data"

# =========================================================
# Create output folders
# =========================================================

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("objects", showWarnings = FALSE)

# =========================================================
# Load datasets
# =========================================================

# Arm_d8
Arm_d8.data <- Read10X_h5(
  file.path(
    data_dir,
    "GSM7919089_SC_Batch_1_Arm_d8_lane_A_sample_feature_bc_matrix.h5"
  )
)

Arm_d8 <- CreateSeuratObject(
  counts = Arm_d8.data$`Gene Expression`,
  project = "Arm_d8",
  min.cells = 3,
  min.features = 200
)

Arm_d8$sample <- "Arm_d8"
Arm_d8$condition <- "Arm"
Arm_d8$day <- 8

# Arm_d21
Arm_d21.data <- Read10X_h5(
  file.path(
    data_dir,
    "GSM7061752_SC_Batch_1_Arm_d21_lane_Asample_feature_bc_matrix.h5"
  )
)

Arm_d21 <- CreateSeuratObject(
  counts = Arm_d21.data$`Gene Expression`,
  project = "Arm_d21",
  min.cells = 3,
  min.features = 200
)

Arm_d21$sample <- "Arm_d21"
Arm_d21$condition <- "Arm"
Arm_d21$day <- 21

# Clone13_d21
Clone13_d21.data <- Read10X_h5(
  file.path(
    data_dir,
    "GSM7061754_SC_Batch_1_Clone13_d21_lane_Asample_feature_bc_matrix.h5"
  )
)

Clone13_d21 <- CreateSeuratObject(
  counts = Clone13_d21.data$`Gene Expression`,
  project = "Clone13_d21",
  min.cells = 3,
  min.features = 200
)

Clone13_d21$sample <- "Clone13_d21"
Clone13_d21$condition <- "Clone13"
Clone13_d21$day <- 21

# =========================================================
# Quality control
# =========================================================

Arm_d8[["percent.mt"]] <- PercentageFeatureSet(
  Arm_d8,
  pattern = "^mt-"
)

Arm_d21[["percent.mt"]] <- PercentageFeatureSet(
  Arm_d21,
  pattern = "^mt-"
)

Clone13_d21[["percent.mt"]] <- PercentageFeatureSet(
  Clone13_d21,
  pattern = "^mt-"
)

# =========================================================
# QC violin plots
# =========================================================

p1 <- VlnPlot(
  Arm_d8,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
) + ggtitle("Arm_d8")

p2 <- VlnPlot(
  Arm_d21,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
) + ggtitle("Arm_d21")

p3 <- VlnPlot(
  Clone13_d21,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
) + ggtitle("Clone13_d21")

(p1 | p2 | p3)

ggsave(
  "figures/QC_violin_plots.png",
  width = 15,
  height = 5,
  dpi = 300
)

# =========================================================
# Filtering
# =========================================================

Arm_d8 <- subset(
  Arm_d8,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 5500 &
    percent.mt < 15
)

Arm_d21 <- subset(
  Arm_d21,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 5500 &
    percent.mt < 15
)

Clone13_d21 <- subset(
  Clone13_d21,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 5500 &
    percent.mt < 15
)

# =========================================================
# Normalization and preprocessing
# =========================================================

objs <- list(
  Arm_d8 = Arm_d8,
  Arm_d21 = Arm_d21,
  Clone13_d21 = Clone13_d21
)

objs <- lapply(objs, function(obj) {
  
  DefaultAssay(obj) <- "RNA"
  
  obj <- NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  
  obj <- ScaleData(
    obj,
    features = VariableFeatures(obj),
    vars.to.regress = "percent.mt",
    verbose = FALSE
  )
  
  obj <- RunPCA(
    obj,
    features = VariableFeatures(obj),
    npcs = 50,
    verbose = FALSE
  )
  
  return(obj)
})

Arm_d8 <- objs$Arm_d8
Arm_d21 <- objs$Arm_d21
Clone13_d21 <- objs$Clone13_d21

# =========================================================
# PCA visualization
# =========================================================

pca1 <- ElbowPlot(Arm_d8, ndims = 50) +
  ggtitle("Arm_d8")

pca2 <- ElbowPlot(Arm_d21, ndims = 50) +
  ggtitle("Arm_d21")

pca3 <- ElbowPlot(Clone13_d21, ndims = 50) +
  ggtitle("Clone13_d21")

(pca1 | pca2 | pca3)

ggsave(
  "figures/ElbowPlots.png",
  width = 12,
  height = 4,
  dpi = 300
)

# =========================================================
# Clustering and UMAP
# =========================================================

dims_use <- 1:30

# Arm_d8
Arm_d8 <- FindNeighbors(
  Arm_d8,
  dims = dims_use,
  verbose = FALSE
)

Arm_d8 <- FindClusters(
  Arm_d8,
  resolution = 0.4,
  verbose = FALSE
)

Arm_d8 <- RunUMAP(
  Arm_d8,
  dims = dims_use,
  min.dist = 0.3,
  verbose = FALSE
)

# Arm_d21
Arm_d21 <- FindNeighbors(
  Arm_d21,
  dims = dims_use,
  verbose = FALSE
)

Arm_d21 <- FindClusters(
  Arm_d21,
  resolution = 0.4,
  verbose = FALSE
)

Arm_d21 <- RunUMAP(
  Arm_d21,
  dims = dims_use,
  min.dist = 0.3,
  verbose = FALSE
)

# Clone13_d21
Clone13_d21 <- FindNeighbors(
  Clone13_d21,
  dims = dims_use,
  verbose = FALSE
)

Clone13_d21 <- FindClusters(
  Clone13_d21,
  resolution = 0.4,
  verbose = FALSE
)

Clone13_d21 <- RunUMAP(
  Clone13_d21,
  dims = dims_use,
  min.dist = 0.3,
  verbose = FALSE
)

# =========================================================
# UMAP plots
# =========================================================

u1 <- DimPlot(
  Arm_d8,
  reduction = "umap",
  label = TRUE
) + ggtitle("Arm_d8")

u2 <- DimPlot(
  Arm_d21,
  reduction = "umap",
  label = TRUE
) + ggtitle("Arm_d21")

u3 <- DimPlot(
  Clone13_d21,
  reduction = "umap",
  label = TRUE
) + ggtitle("Clone13_d21")

(u1 | u2 | u3)

ggsave(
  "figures/UMAP_three_samples.png",
  width = 15,
  height = 5,
  dpi = 300
)

# =========================================================
# Merge datasets
# =========================================================

combined <- merge(
  Arm_d8,
  y = c(Arm_d21, Clone13_d21),
  add.cell.ids = c(
    "Arm_d8",
    "Arm_d21",
    "Clone13_d21"
  )
)

DefaultAssay(combined) <- "RNA"

combined <- NormalizeData(
  combined,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

combined <- FindVariableFeatures(
  combined,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

combined <- ScaleData(
  combined,
  features = VariableFeatures(combined),
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

combined <- RunPCA(
  combined,
  features = VariableFeatures(combined),
  npcs = 50,
  verbose = FALSE
)

combined <- FindNeighbors(
  combined,
  dims = dims_use,
  verbose = FALSE
)

combined <- FindClusters(
  combined,
  resolution = 0.4,
  verbose = FALSE
)

combined <- RunUMAP(
  combined,
  dims = dims_use,
  min.dist = 0.3,
  verbose = FALSE
)

# =========================================================
# IMPORTANT Seurat v5 fix
# Join layers BEFORE marker analysis
# =========================================================

combined <- JoinLayers(combined)

# =========================================================
# Combined UMAPs
# =========================================================

p_sample <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "orig.ident"
) + ggtitle("Combined UMAP by Sample")

p_cluster <- DimPlot(
  combined,
  reduction = "umap",
  label = TRUE
) + ggtitle("Combined UMAP by Cluster")

(p_sample | p_cluster)

ggsave(
  "figures/Combined_UMAP.png",
  width = 12,
  height = 5,
  dpi = 300
)

# =========================================================
# Marker analysis
# =========================================================

Idents(combined) <- combined$seurat_clusters

combined_markers <- FindAllMarkers(
  combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

write.csv(
  combined_markers,
  "results/combined_cluster_markers.csv",
  row.names = FALSE
)

# =========================================================
# Top markers per cluster
# =========================================================

fc_col <- if ("avg_log2FC" %in% colnames(combined_markers)) {
  "avg_log2FC"
} else {
  "avg_logFC"
}

cluster_col <- if ("cluster" %in% colnames(combined_markers)) {
  "cluster"
} else {
  "ident"
}

top10 <- combined_markers %>%
  group_by(.data[[cluster_col]]) %>%
  slice_max(order_by = .data[[fc_col]], n = 10)

write.csv(
  top10,
  "results/combined_cluster_markers_top10.csv",
  row.names = FALSE
)

# =========================================================
# Marker visualization
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

dotplot <- DotPlot(
  combined,
  features = genes
) + RotatedAxis()

print(dotplot)

ggsave(
  "figures/Tcell_marker_dotplot.png",
  width = 12,
  height = 6,
  dpi = 300
)

# =========================================================
# Manual cluster annotation
# =========================================================

map <- c(
  "3" = "Memory T cells",
  "2" = "Effector T cells",
  "5" = "Tpex cells",
  "4" = "Effector T cells",
  "1" = "Exhausted T cells"
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
# Save final object
# =========================================================

saveRDS(
  combined,
  file = "objects/combined_integrated_annotated_final.rds"
)

message("Analysis complete.")