# =============================================================================
# FIGURE 1: Single-Cell RNA-Seq Analysis of GSE115978 Melanoma
# Panels: A (UMAP by cluster), B (clustree), C (UMAP by sample), D (top markers heatmap)
# =============================================================================

# 0. SETUP
# -----------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "Seurat", "dplyr", "scater", "pheatmap", "scales",
  "clustree", "ComplexHeatmap", "circlize")

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs, ask = FALSE)
lapply(required_pkgs, library, character.only = TRUE)

options(stringsAsFactors = FALSE)
theme_set(theme_classic())

unicode_minus <- function(x) sub("^-", "\\u2212", x)

# 1. DATA LOADING (MALIGNANT CELLS ONLY)
# -----------------------------------------------------------------------------
cat("Loading GSE115978 data...\n")

counts <- read.csv("GSE115978_counts.csv", header = TRUE, sep = ",", row.names = 1)
annotation <- read.csv("GSE115978_annotations.csv", header = TRUE, sep = ",")

# Malignant cells
annotation_mal <- annotation[annotation$cell.types == "Malignant", ]

primary_samples <- c("Mel84", "Mel129pa", "Mel129pb", "Mel105")
metastatic_samples <- c("Mel53", "Mel71", "Mel79", "Mel80", "Mel81", "Mel82",
                        "Mel103", "Mel116", "Mel112", "Mel478", "Mel128")

primary_cells <- annotation_mal$cells[annotation_mal$samples %in% primary_samples]
metastatic_cells <- annotation_mal$cells[annotation_mal$samples %in% metastatic_samples]

counts_primary <- counts[, primary_cells]
counts_metastatic <- counts[, metastatic_cells]

anno_primary <- annotation_mal[annotation_mal$cells %in% primary_cells, ]
anno_metastatic <- annotation_mal[annotation_mal$cells %in% metastatic_cells, ]

# 2. SEURAT OBJECT CREATION AND QC
# -----------------------------------------------------------------------------
cat("Creating Seurat objects...\n")

seurat_primary <- CreateSeuratObject(
  counts = counts_primary,
  project = "Primary",
  meta.data = data.frame(
    cells    = anno_primary$cells,
    patient  = anno_primary$samples,
    type     = anno_primary$type,
    celltype = anno_primary$cell.types,
    row.names = anno_primary$cells
  )
)

seurat_metastatic <- CreateSeuratObject(
  counts = counts_metastatic,
  project = "Metastatic",
  meta.data = data.frame(
    cells    = anno_metastatic$cells,
    patient  = anno_metastatic$samples,
    type     = anno_metastatic$type,
    celltype = anno_metastatic$cell.types,
    row.names = anno_metastatic$cells
  )
)


# Merge objects
seurat <- merge(
  seurat_primary,
  seurat_metastatic,
  add.cell.ids = c("Primary", "Metastatic"),
  project = "GSE115978"
)

# Calculate QC metrics
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, pattern = "^RP[SL]")

# Filter high-quality cells (nFeature_RNA > 2500)
seurat <- subset(seurat, subset = nFeature_RNA > 2500)


# 3. NORMALIZATION, PCA, UMAP
# -----------------------------------------------------------------------------
cat("Normalization, PCA, UMAP...\n")

seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(seurat, features = rownames(seurat))

seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
seurat <- FindNeighbors(seurat, dims = 1:5)
seurat <- RunUMAP(seurat, dims = 1:5)

# 4. RESOLUTION OPTIMIZATION WITH CLUSTREE (FIG 1B)
# -----------------------------------------------------------------------------
cat("Optimizing resolution with clustree...\n")

resolutions <- c(0, 0.2, 0.4, 0.8, 1.2)
seurat_res <- seurat

# Run clustering at multiple resolutions
seurat_res <- FindClusters(seurat_res, resolution = resolutions)

# Add UMAP coordinates to metadata for clustree_overlay
umap_coords <- as.data.frame(Embeddings(seurat_res[["umap"]]))
metadata <- seurat_res@meta.data
metadata$UMI_id <- rownames(metadata)
umap_coords$UMI_id <- rownames(umap_coords)
metadata_merged <- merge(metadata, umap_coords, by = "UMI_id")
rownames(metadata_merged) <- metadata_merged$UMI_id
metadata_merged$UMI_id <- NULL
seurat_res@meta.data <- metadata_merged

# Panel B: clustree
p_tree <- clustree(seurat_res, node_size = 12, node_text_size = 8) +
  theme(
    legend.text  = element_text(size = 13, colour = "black", face = "bold"),
    legend.title = element_text(size = 13, colour = "black", face = "bold")
  )

p_umap_res <- clustree_overlay(
  seurat_res,
  x_value = "UMAP_1",
  y_value = "UMAP_2",
  red_dim = "umap",
  suffix = "UMAP"
)

ggsave("Fig1B_clustree_tree.png", p_tree, width = 8, height = 6.5, dpi = 1200)
ggsave("Fig1B_clustree_umap_overlay.png", p_umap_res, width = 8, height = 7, dpi = 1200)

cat("Selected resolution = 0.4 (6 stable, biologically interpretable clusters).\n")

# 5. FINAL CLUSTERING AT RESOLUTION 0.4
# -----------------------------------------------------------------------------
seurat <- FindClusters(seurat, resolution = 0.4)
Idents(seurat) <- "seurat_clusters"

# Order cluster levels (0,5,1,4,3,2) to match your original script/figure
levels(seurat) <- c("0", "5", "1", "4", "3", "2")

# 6. UMAP PANELS (FIG 1A AND 1C)
# -----------------------------------------------------------------------------

## Panel A: UMAP by cluster (0–5)
cat("Generating UMAP by cluster (Fig 1A)...\n")

p_umap_cluster <- DimPlot(
  seurat,
  reduction = "umap",
  group.by  = "seurat_clusters",
  label     = TRUE,
  label.size = 6,
  repel     = TRUE,
  pt.size   = 1.0
) +
  theme(
    plot.title  = element_blank(),
    axis.title  = element_text(size = 20, colour = "black", face = "bold"),
    axis.text   = element_text(size = 14, colour = "black", face = "bold"),
    legend.text = element_text(size = 12, colour = "black", face = "bold")
  ) +
  scale_x_continuous(labels = unicode_minus) +
  scale_y_continuous(labels = unicode_minus)

ggsave("Fig1A_UMAP_clusters.png", p_umap_cluster, width = 6, height = 6, dpi = 1200)

## Panel C: UMAP by patient/sample
cat("Generating UMAP by sample (Fig 1C)...\n")

p_umap_patient <- DimPlot(
  seurat,
  reduction = "umap",
  group.by  = "patient",
  label     = FALSE,
  pt.size   = 1.0
) +
  theme(
    plot.title  = element_blank(),
    axis.title  = element_text(size = 20, colour = "black", face = "bold"),
    axis.text   = element_text(size = 14, colour = "black", face = "bold"),
    legend.text = element_text(size = 10, colour = "black", face = "bold")
  ) +
  scale_x_continuous(labels = unicode_minus) +
  scale_y_continuous(labels = unicode_minus)

ggsave("Fig1C_UMAP_patient.png", p_umap_patient, width = 6, height = 6, dpi = 1200)

# 7. DIFFERENTIAL EXPRESSION AND HEATMAP (FIG 1D)
# -----------------------------------------------------------------------------
cat("Differential expression and heatmap (Fig 1D)...\n")

deg_markers <- FindAllMarkers(
  seurat,
  only.pos       = TRUE,
  min.pct        = 0.25,
  logfc.threshold = 1
)

top10_markers <- deg_markers %>%
  group_by(cluster) %>%
  top_n(10, wt = avg_log2FC)

# Heatmap of top 10 markers per cluster
p_heatmap <- DoHeatmap(
  seurat,
  features = top10_markers$gene,
  angle    = 0,
  label    = FALSE
) +
  theme(
    axis.text.y   = element_text(size = 8, colour = "black", face = "bold.italic"),
    legend.key.size = unit(1.0, "cm"),
    legend.text   = element_text(size = 10, face = "bold"),
    legend.title  = element_text(size = 10, face = "bold")
  ) +
  scale_color_discrete(name = "Cluster", na.translate = FALSE)

ggsave("Fig1D_top10_markers_heatmap.png", p_heatmap, width = 8, height = 10, dpi = 1200)

# Save marker tables
write.csv(top10_markers, "Fig1_top10_markers_per_cluster.csv", row.names = FALSE)
write.csv(deg_markers, "Fig1_all_DEGs_logFC1.csv", row.names = FALSE)

cat("FIGURE 1 SCRIPT FINISHED.\n")
