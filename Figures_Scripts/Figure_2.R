# =============================================================================
# FIGURE 2: Pathway Enrichment of GSE115978 Melanoma Clusters
# Panels:
#   A – GO BP enrichment (clusterProfiler::compareCluster dotplot)
#   B – Hallmark GSVA heatmap (pheatmap)
#   C – IFN pathway gene expression heatmap (ComplexHeatmap)
# Inputs:
#   - Seurat object with malignant cells and clusters at resolution 0.4
#   - DEG.Markers (FindAllMarkers result, logFC >= 1)
# =============================================================================

# 0. SETUP ---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "Seurat", "dplyr", "clusterProfiler", "org.Hs.eg.db", "scales",
  "escape", "dittoSeq", "UCell", "GSVA",
  "pheatmap", "grid", "gridExtra", "circlize", "ComplexHeatmap"
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs, ask = FALSE)
lapply(required_pkgs, library, character.only = TRUE)

options(stringsAsFactors = FALSE)

# Assumes objects already exist in environment:
#   seurat        – Seurat object from Figure 1 (resolution 0.4, malignant cells only)
#   DEG.Markers   – output of FindAllMarkers(seurat, only.pos=TRUE, logfc.threshold=1)
#   (In your original script, this was GSE115978.prim.and.metast)


# 1. GO ENRICHMENT (FIG 2A) ----------------------------------------------------
cat("Running GO BP enrichment (Fig 2A)...\n")

# Split DEG markers by cluster
cluster0.genes <- DEG.Markers$gene[DEG.Markers$cluster == 0]
cluster5.genes <- DEG.Markers$gene[DEG.Markers$cluster == 5]
cluster1.genes <- DEG.Markers$gene[DEG.Markers$cluster == 1]
cluster4.genes <- DEG.Markers$gene[DEG.Markers$cluster == 4]
cluster3.genes <- DEG.Markers$gene[DEG.Markers$cluster == 3]
cluster2.genes <- DEG.Markers$gene[DEG.Markers$cluster == 2]

list.for.GO <- list(
  "0" = cluster0.genes,
  "5" = cluster5.genes,
  "1" = cluster1.genes,
  "4" = cluster4.genes,
  "3" = cluster3.genes,
  "2" = cluster2.genes
)

GO.res <- compareCluster(
  geneCluster = list.for.GO,
  fun        = enrichGO,
  OrgDb      = "org.Hs.eg.db",
  keyType    = "SYMBOL",
  ont        = "BP"
)

write.csv(GO.res, "Fig2A_GO_Enrichment_logFC1.csv", row.names = FALSE)

png("Fig2A_GO_dotplot.png", width = 11, height = 9, units = "in", res = 1200)
clusterProfiler::dotplot(GO.res) +
  theme(
    plot.title  = element_text(size = 22, colour = "black", face = "bold", hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 16, colour = "black", face = "bold"),
    axis.text.x = element_text(size = 25, colour = "black", face = "bold"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 17, face = "bold")
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  scale_x_discrete(labels = c("0", "5", "1", "4", "3", "2"))
dev.off()


# 2. HALLMARK GSVA HEATMAP (FIG 2B) -------------------------------------------
cat("Running GSVA Hallmark enrichment (Fig 2B)...\n")

gene.sets.H <- getGeneSets(library = "H", species = "Homo sapiens")

# Average expression per cluster
Average.GSE115978 <- AverageExpression(GSE115978.prim.and.metast, return.seurat = TRUE)
avg_mat <- Average.GSE115978[["RNA"]]@data

# (Optional) restrict to DEG genes only
avg_mat_sub <- avg_mat[DEG.Markers$gene, , drop = FALSE]

# Full matrix is used for GSVA (matches your original code)
GSVA.Enrich <- gsva(
  avg_mat,
  gene.sets.H,
  min.sz = 5,
  max.sz = 500,
  kcdf  = "Gaussian"
)

write.csv(GSVA.Enrich, "Fig2B_Hallmark_GSVA_scores.csv")

ROWnames <- lapply(
  substr(gsub("_", " ", gsub("", "", rownames(GSVA.Enrich))), 1, 500),
  function(x) bquote(bold(.(x)))
)
COLnames <- lapply(
  colnames(GSVA.Enrich),
  function(x) bquote(bold(.(x)))
)

png("Fig2B_HALLMARKS_heatmap.png", width = 10, height = 9, units = "in", res = 1200)
pheatmap::pheatmap(
  GSVA.Enrich,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  angle_col         = 0,
  scale             = "column",
  clustering_method = "ward.D2",
  fontsize          = 14,
  labels_col        = as.expression(COLnames),
  labels_row        = as.expression(ROWnames)
)
dev.off()

# 3. IFN PATHWAY HEATMAP (FIG 2C) ---------------------------------------------
cat("Generating IFN pathway heatmap (Fig 2C)...\n")

# Hallmark IFN gene sets
Hallmark.IFN.alpha  <- gene.sets.H[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]]@geneIds
Hallmark.IFN.gamma  <- gene.sets.H[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]@geneIds

# Use scaled average expression per cluster
Scaled.count <- as.matrix(Average.GSE115978[["RNA"]]@scale.data)
# Select clusters 5, 1, 3 (column order: 5, 1, 3)
Scaled.count <- Scaled.count[, c("5", "1", "3")]

Scaled.count.IFNa.5.1.3 <- Scaled.count[rownames(Scaled.count) %in% Hallmark.IFN.alpha, , drop = FALSE]
Scaled.count.IFNg.5.1.3 <- Scaled.count[rownames(Scaled.count) %in% Hallmark.IFN.gamma, , drop = FALSE]

# Color scheme
my.Color  <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
my.Breaks <- seq(-3, 3, length.out = 100)

IFN.heatmap <- ComplexHeatmap::Heatmap(
  Scaled.count.IFNg.5.1.3,
  column_names_rot   = 360,
  column_title       = NULL,
  name               = "Expression",
  show_column_names  = TRUE,
  cluster_columns    = FALSE,
  show_column_dend   = FALSE,
  show_row_dend      = FALSE,
  show_row_names     = TRUE,
  col                = circlize::colorRamp2(my.Breaks, my.Color),
  row_names_gp       = grid::gpar(fontsize = 12, fontface = "bold.italic"),
  column_names_gp    = grid::gpar(fontsize = 18, fontface = "bold"),
  heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 14, fontface = "bold"))
)

png("Fig2C_IFNg_heatmap.png", width = 6, height = 12, units = "in", res = 1200)
draw(IFN.heatmap)
dev.off()

cat("FIGURE 2 SCRIPT FINISHED.\n")
cat("Generated files:\n")
cat("- Fig2A_GO_dotplot.png\n")
cat("- Fig2B_HALLMARKS_heatmap.png\n")
cat("- Fig2C_IFNg_heatmap.png\n")
cat("- Fig2A_GO_Enrichment_logFC1.csv\n")
cat("- Fig2B_Hallmark_GSVA_scores.csv\n")
