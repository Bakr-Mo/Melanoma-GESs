# =============================================================================
# FIGURE 3: Immune–Melanoma Cell–Cell Communication (CellChat)
# Panels:
#   A – UMAP of immune cell types (Monaco main annotation)
#   B – Aggregated cell–cell communication circle plot (immune + melanoma)
#   C – Ligand–receptor bubble plot for selected interactions
#   D/E – MHC-I / MHC-II pathway networks (signaling to T cells / melanoma)
#   F – MHC-related chord diagram
#   G – IFN-II signaling pathway network / roles
# =============================================================================

# 0. SETUP ---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "Seurat", "dplyr", "ggplot2", "SingleR", "celldex",
  "CellChat", "patchwork", "NMF", "ggalluvial",
  "future"
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs, ask = FALSE)
lapply(required_pkgs, library, character.only = TRUE)

options(stringsAsFactors = FALSE)
theme_set(theme_classic())
set.seed(1234)

# Load malignant object which have been already saved from Figure 1, load it:
# GSE115978.malignant <- readRDS("GSE115978_malignant_seurat.rds")

# 1. LOAD DATA AND SUBSET IMMUNE CELLS ----------------------------------------
cat("Loading GSE115978 and subsetting immune cells...\n")

counts_all <- read.csv("GSE115978_counts.csv", header = TRUE, sep = ",", row.names = 1)
annot_all  <- read.csv("GSE115978_annotations.csv", header = TRUE, sep = ",")

immune_types <- c("T.cell", "T.CD4", "T.CD8", "B.cell", "Macrophage", "NK")
immune_samples <- c(
  "Mel84", "Mel129pa", "Mel129pb", "Mel105",
  "Mel53", "Mel71", "Mel79", "Mel80", "Mel81", "Mel82",
  "Mel103", "Mel116", "Mel112", "Mel478", "Mel128"
)

annot_immune <- annot_all[
  annot_all$samples %in% immune_samples &
    annot_all$cell.types %in% immune_types,
]

counts_immune <- counts_all[, annot_immune$cells]

immune_seurat <- CreateSeuratObject(
  counts = counts_immune,
  project = "GSE115978_Immune",
  meta.data = data.frame(
    cells    = annot_immune$cells,
    patient  = annot_immune$samples,
    type     = annot_immune$type,
    celltype = annot_immune$cell.types,
    row.names = annot_immune$cells
  )
)

# 2. SEURAT PIPELINE FOR IMMUNE CELLS -----------------------------------------
cat("Running Seurat pipeline for immune cells...\n")

immune_seurat[["percent.mt"]] <- PercentageFeatureSet(immune_seurat, pattern = "^MT-")
immune_seurat[["percent.rb"]] <- PercentageFeatureSet(immune_seurat, pattern = "^RP[SL]")

immune_seurat <- subset(immune_seurat, subset = nFeature_RNA > 2500)

immune_seurat <- NormalizeData(immune_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
immune_seurat <- FindVariableFeatures(immune_seurat, selection.method = "vst", nfeatures = 2000)
immune_seurat <- ScaleData(immune_seurat, features = rownames(immune_seurat))

immune_seurat <- RunPCA(immune_seurat, features = VariableFeatures(immune_seurat))
immune_seurat <- FindNeighbors(immune_seurat, dims = 1:5)
immune_seurat <- FindClusters(immune_seurat, resolution = 0.8)
immune_seurat <- RunUMAP(immune_seurat, dims = 1:5)

# 3. IMMUNE CELL ANNOTATION WITH MONACO (FIG 3A) ------------------------------
cat("Annotating immune cells with Monaco (SingleR)...\n")

monaco_ref <- celldex::MonacoImmuneData()
sce_immune <- as.SingleCellExperiment(DietSeurat(immune_seurat))

monaco_main <- SingleR(
  test          = sce_immune,
  assay.type.test = 1,
  ref           = monaco_ref,
  labels        = monaco_ref$label.main
)

immune_seurat$monaco.main <- monaco_main$pruned.labels

# Replace remaining NA with CD4+ T cells (as in your original script)
immune_seurat$monaco.main[is.na(immune_seurat$monaco.main)] <- "CD4+ T cells"

Idents(immune_seurat) <- "monaco.main"

p3A <- DimPlot(
  immune_seurat,
  reduction = "umap",
  group.by  = "monaco.main",
  label     = TRUE,
  repel     = TRUE,
  label.size = 3,
  pt.size   = 1.2
) +
  theme(
    plot.title  = element_blank(),
    axis.title  = element_text(size = 15, face = "bold", colour = "black"),
    axis.text   = element_text(size = 12, face = "bold", colour = "black"),
    legend.text = element_text(size = 11, face = "bold", colour = "black")
  )

ggsave("Fig3A_UMAP_Immune_MonacoMain.png", p3A, width = 6, height = 5, dpi = 700)

# 4. MERGE MALIGNANT + IMMUNE AND PREPARE FOR CELLCHAT ------------------------
cat("Merging malignant + immune for CellChat...\n")

# Ensure malignant object exists. If not, stop here.
if (!exists("GSE115978.malignant"))
  stop("GSE115978.malignant Seurat object not found. Load from Figure 1 output.")

GSE115978.malignant$ident <- Idents(GSE115978.malignant)
immune_seurat$ident        <- Idents(immune_seurat)

malig_immune <- merge(
  GSE115978.malignant,
  immune_seurat,
  merge.data   = TRUE,
  add.cell.ids = c("Malignant", "Immune"),
  project      = "GSE115978_Malig_Immune"
)

Idents(malig_immune) <- "ident"

# 5. CELlCHAT OBJECT AND CORE INFERENCE ----------------------------------------
cat("Creating CellChat object and running core inference...\n")

cellchat.Immune <- createCellChat(
  object = malig_immune,
  meta   = malig_immune@meta.data,
  group.by = "ident"
)

CellChatDB      <- CellChatDB.human
CellChatDB.use  <- CellChatDB
cellchat.Immune@DB <- CellChatDB.use

cellchat.Immune <- subsetData(cellchat.Immune)

future::plan("multisession", workers = 6)

cellchat.Immune <- identifyOverExpressedGenes(cellchat.Immune)
cellchat.Immune <- identifyOverExpressedInteractions(cellchat.Immune)
cellchat.Immune <- projectData(cellchat.Immune, PPI.human)

cellchat.Immune <- computeCommunProb(cellchat.Immune)
cellchat.Immune <- computeCommunProbPathway(cellchat.Immune)
cellchat.Immune <- aggregateNet(cellchat.Immune)

group_size <- as.numeric(table(cellchat.Immune@idents))

# 6. FIGURE 3B – AGGREGATED NETWORK CIRCLE PLOT -------------------------------
cat("Plotting aggregated communication circle plot (Fig 3B)...\n")

png("Fig3B_CellChat_Aggregated_Count.png", width = 6, height = 6, units = "in", res = 700)
netVisual_circle(
  cellchat.Immune@net$count,
  vertex.weight = group_size,
  weight.scale  = TRUE,
  label.edge    = FALSE,
  title.name    = "Number of interactions",
  margin        = c(0.2, 0.2, 0.2, 0.2),
  vertex.label.cex = 0.9
)
dev.off()

# 7. FIGURE 3C – LIGAND–RECEPTOR BUBBLE PLOT ----------------------------------
cat("Plotting ligand–receptor bubble plot (Fig 3C)...\n")

# Adjust sources/targets indices to match your idents order
levels(cellchat.Immune@idents)

# Example: malignant as source to CD4+/CD8+ T cells as targets
# Replace indices below with your actual positions
sources_use <- which(levels(cellchat.Immune@idents) %in% c("T cells", "Melanogenesis", "EMT"))[1]
targets_use <- which(levels(cellchat.Immune@idents) %in% c("CD4+ T cells", "CD8+ T cells", "T cells"))

png("Fig3C_LR_Bubble_Malig_to_Tcells.png", width = 4, height = 10, units = "in", res = 700)
netVisual_bubble(
  cellchat.Immune,
  sources.use   = sources_use,
  targets.use   = targets_use,
  remove.isolate = FALSE,
  thresh       = 0.001
)
dev.off()

# 8. MHC-I / MHC-II / IFN-II NETWORKS (FIG 3D–G STYLE) ------------------------
cat("Plotting MHC-I / MHC-II / IFN-II signaling networks...\n")

# Check available pathways
pathways_all <- cellchat.Immune@netP$pathways

# Example pathways: "MHC-I", "MHC-II", "IFN-II"
pathway_MHCI   <- "MHC-I"
pathway_MHCII  <- "MHC-II"
pathway_IFNII  <- "IFN-II"

# Choose T-cell receiving vertices (update indices based on levels)
vertex_receiver_T <- which(levels(cellchat.Immune@idents) %in% c("CD4+ T cells", "CD8+ T cells", "T cells"))

## MHC-I hierarchy plot (similar to 3D/3E)
png("Fig3D_MHCI_Hierarchy_Tcells_vs_Melanoma.png", width = 7, height = 5, units = "in", res = 700)
netVisual_aggregate(
  cellchat.Immune,
  signaling       = pathway_MHCI,
  vertex.receiver = vertex_receiver_T,
  layout          = "hierarchy",
  thresh          = 0.01
)
dev.off()

## MHC-II hierarchy plot
png("Fig3E_MHCII_Hierarchy_Tcells_vs_Melanoma.png", width = 7, height = 5, units = "in", res = 700)
netVisual_aggregate(
  cellchat.Immune,
  signaling       = pathway_MHCII,
  vertex.receiver = vertex_receiver_T,
  layout          = "hierarchy",
  thresh          = 0.01
)
dev.off()

# 9. CHORD DIAGRAM FOR MHC-I/II (FIG 3F-LIKE) ---------------------------------
cat("Plotting chord diagram for MHC-I/II...\n")

png("Fig3F_MHC_I_II_Chord_Tcells.png", width = 7, height = 7, units = "in", res = 700)
par(mar = c(0, 0, 0, 0))
netVisual_chord_gene(
  cellchat.Immune,
  sources.use = vertex_receiver_T,
  targets.use = vertex_receiver_T,
  signaling   = c(pathway_MHCI, pathway_MHCII),
  legend.pos.x = 10,
  legend.pos.y = 10,
  small.gap   = 4,
  big.gap     = 5
)
dev.off()

# 10. IFN-II SIGNALLING ROLE / NETWORK (FIG 3G STYLE) -------------------------
cat("Plotting IFN-II signaling roles (Fig 3G)...\n")

cellchat.Immune <- netAnalysis_computeCentrality(cellchat.Immune, slot.name = "netP")

png("Fig3G_IFNII_SignalingRole_Network.png", width = 5, height = 4, units = "in", res = 700)
netAnalysis_signalingRole_network(
  cellchat.Immune,
  signaling = pathway_IFNII,
  width     = 8,
  height    = 3,
  font.size = 10
)
dev.off()

# Optional: scatter plot for roles
gg_all <- netAnalysis_signalingRole_scatter(cellchat.Immune)
gg_ifn <- netAnalysis_signalingRole_scatter(cellchat.Immune, signaling = pathway_IFNII)

ggsave("Fig3G_IFNII_SignalingRole_Scatter_All.png", gg_all + gg_ifn, width = 10, height = 5, dpi = 700)

# 11. SAVE OBJECT -------------------------------------------------------------
saveRDS(cellchat.Immune, file = "cellchat_MaligImmune_Fig3.rds")

cat("FIGURE 3 SCRIPT FINISHED.\n")
cat("Generated key panels:\n")
cat("- Fig3A_UMAP_Immune_MonacoMain.png\n")
cat("- Fig3B_CellChat_Aggregated_Count.png\n")
cat("- Fig3C_LR_Bubble_Malig_to_Tcells.png\n")
cat("- Fig3D_MHCI_Hierarchy_Tcells_vs_Melanoma.png\n")
cat("- Fig3E_MHCII_Hierarchy_Tcells_vs_Melanoma.png\n")
cat("- Fig3F_MHC_I_II_Chord_Tcells.png\n")
cat("- Fig3G_IFNII_SignalingRole_*.png\n")

