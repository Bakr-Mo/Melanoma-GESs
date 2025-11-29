# =============================================================================
# FIGURE 4: TCGA-SKCM Classification by Melanoma GESs
# Panels:
#   A – Stacked bar chart of GES fractions per TCGA-SKCM patient (CIBERSORTx)
#   B – Pie chart of overall GES proportions
#   C – Hallmark GSVA heatmap across TCGA GESs
#   D – xCell cell-type enrichment heatmap across TCGA GESs
# =============================================================================

# 0. SETUP ---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "TCGAbiolinks", "SummarizedExperiment", "ggplot2", "limma", "biomaRt",
  "tibble", "edgeR", "dplyr", "tidyr", "scales", "data.table",
  "Seurat", "msigdbr", "fgsea", "escape", "dittoSeq", "UCell", "GSVA",
  "pheatmap", "xCell"
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs, ask = FALSE)
lapply(required_pkgs, library, character.only = TRUE)

options(stringsAsFactors = FALSE)
theme_set(theme_classic())
set.seed(1234)

# 1. DOWNLOAD AND PREPARE TCGA-SKCM COUNTS ------------------------------------
cat("Downloading TCGA-SKCM HTSeq-Counts with TCGAbiolinks...\n")

query.htseq <- GDCquery(
  project      = "TCGA-SKCM",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts"
)

GDCdownload(query.htseq)
data <- GDCprepare(query.htseq)  # SummarizedExperiment

save(data, file = "SKCMhg38.rda")  # optional

counts_tcga <- assay(data)
clinical    <- data@colData

# 2. ANNOTATE GENES WITH BIOMART AND GET PROTEIN-CODING logCPM ----------------
cat("Annotating genes with biomaRt and computing logCPM...\n")

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "m.ensembl.org")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

attributeNames <- c("ensembl_gene_id", "external_gene_name", "gene_biotype", "description")
annot_all <- getBM(
  attributes = attributeNames,
  filters    = "ensembl_gene_id",
  values     = rownames(counts_tcga),
  mart       = ensembl
)

annot_all_genes <- as.data.frame(counts_tcga) %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(annot_all, by = "ensembl_gene_id")

all.genes.pc <- annot_all_genes[annot_all_genes$gene_biotype == "protein_coding", ]
all.genes.pc <- all.genes.pc %>% distinct(external_gene_name, .keep_all = TRUE)

all_gene_names <- all.genes.pc$external_gene_name
count_matrix_pc <- as.matrix(all.genes.pc[, 2:(ncol(all.genes.pc) - 3)])

logCPM.all <- edgeR::cpm(count_matrix_pc, prior.count = 2, log = TRUE)
rownames(logCPM.all) <- all_gene_names

fwrite(as.data.frame(logCPM.all), "TCGA_SKCM_logCPM.txt", sep = "\t", row.names = TRUE)

# 3. LOAD CIBERSORTx FRACTIONS AND MAP GES CLUSTERS ---------------------------
cat("Loading CIBERSORTx fractions...\n")
CIBERSORTx <- read.csv("CIBERSORTx.csv", row.names = 1, header = TRUE)

CIBERSORTxFractions <- CIBERSORTx[, 1:6]
colnames(CIBERSORTxFractions) <- c("0", "1", "2", "3", "4", "5")

CIBERSORTxFractions$dominant.cluster <- colnames(CIBERSORTxFractions)[
  max.col(CIBERSORTxFractions, ties.method = "first")
]

CIBERSORTx.TCGA.Cl0 <- rownames(CIBERSORTxFractions[CIBERSORTxFractions$dominant.cluster == "0", ])
CIBERSORTx.TCGA.Cl1 <- rownames(CIBERSORTxFractions[CIBERSORTxFractions$dominant.cluster == "1", ])
CIBERSORTx.TCGA.Cl2 <- rownames(CIBERSORTxFractions[CIBERSORTxFractions$dominant.cluster == "2", ])
CIBERSORTx.TCGA.Cl3 <- rownames(CIBERSORTxFractions[CIBERSORTxFractions$dominant.cluster == "3", ])
CIBERSORTx.TCGA.Cl4 <- rownames(CIBERSORTxFractions[CIBERSORTxFractions$dominant.cluster == "4", ])
CIBERSORTx.TCGA.Cl5 <- rownames(CIBERSORTxFractions[CIBERSORTxFractions$dominant.cluster == "5", ])

clinical.w.CIBERSORTx <- as.data.frame(clinical)
clinical.w.CIBERSORTx$scRNAcluster <- ifelse(
  rownames(clinical.w.CIBERSORTx) %in% CIBERSORTx.TCGA.Cl2, "EMT",
  ifelse(rownames(clinical.w.CIBERSORTx) %in% CIBERSORTx.TCGA.Cl0, "Anti-apoptosis",
         ifelse(rownames(clinical.w.CIBERSORTx) %in% CIBERSORTx.TCGA.Cl1, "Melanogenesis",
                ifelse(rownames(clinical.w.CIBERSORTx) %in% CIBERSORTx.TCGA.Cl4, "Ribosomal biogenesis",
                       ifelse(rownames(clinical.w.CIBERSORTx) %in% CIBERSORTx.TCGA.Cl5, "Immune cell interactions",
                              ifelse(rownames(clinical.w.CIBERSORTx) %in% CIBERSORTx.TCGA.Cl3, "Extracellular structure organization",
                                     "Others"))))))
)

clinical.w.CIBERSORTx$scRNAcluster <- factor(
  clinical.w.CIBERSORTx$scRNAcluster,
  levels = c(
    "Anti-apoptosis",
    "Immune cell interactions",
    "Melanogenesis",
    "Ribosomal biogenesis",
    "Extracellular structure organization",
    "EMT"
  ),
  ordered = TRUE
)

# 4. FIGURE 4A – STACKED BARPLOT (CIBERSORTx FRACTIONS) -----------------------
cat("Generating Fig 4A stacked barplot...\n")

CIBERSORTxFractions_for_plot <- CIBERSORTx[, 1:6]
CIBERSORTxFractions_for_plot <- CIBERSORTxFractions_for_plot[, c(1, 6, 2, 5, 4, 3)]
colnames(CIBERSORTxFractions_for_plot) <- c(
  "0 : Anti-apoptosis",
  "5 : Immune cell interactions",
  "1 : Melanogenesis",
  "4 : Ribosomal biogenesis",
  "3 : Extracellular structure organization",
  "2 : EMT"
)

CIBERSORTxFractions_for_plot$mixture <- rownames(CIBERSORTxFractions_for_plot)
rownames(CIBERSORTxFractions_for_plot) <- NULL

test_df <- CIBERSORTxFractions_for_plot %>%
  pivot_longer(!mixture, names_to = "cell_types", values_to = "percent") %>%
  rename(barcodes = mixture)

test_df$cell_types <- factor(
  test_df$cell_types,
  levels = c(
    "0 : Anti-apoptosis",
    "5 : Immune cell interactions",
    "1 : Melanogenesis",
    "4 : Ribosomal biogenesis",
    "3 : Extracellular structure organization",
    "2 : EMT"
  ),
  ordered = TRUE
)

p4A <- ggplot(test_df, aes(x = barcodes, y = percent, fill = cell_types)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(fill = "Clusters : GESs") +
  geom_col(width = 1) +
  scale_fill_manual(values = c("black", "green", "blue", "yellow", "orange", "brown")) +
  labs(x = "TCGA-SKCM Patients", y = "Relative Percent") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 20, colour = "black", face = "bold"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 20, colour = "black", face = "bold"),
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(size = 13, colour = "black", face = "bold"),
    legend.text  = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 17, face = "bold")
  )

ggsave("Fig4A_TCGA_CIBERSORTx_StackedBar_GESs.png", p4A, width = 12, height = 6, dpi = 1200)

# 5. FIGURE 4B – PIE CHART OF GES PROPORTIONS ---------------------------------
cat("Generating Fig 4B pie chart...\n")

Clusters.Percentage <- data.frame(
  GESs = c(
    "Anti-apoptosis",
    "Immune cell interactions",
    "Melanogenesis",
    "Ribosomal biogenesis",
    "Extracellular structure organization",
    "EMT"
  ),
  Percentage = c(5.74, 32.55, 3.83, 27.45, 27.45, 2.98)
)

Clusters.Percentage$GESs <- factor(
  Clusters.Percentage$GESs,
  ordered = TRUE,
  levels = c(
    "Anti-apoptosis",
    "Immune cell interactions",
    "Melanogenesis",
    "Ribosomal biogenesis",
    "Extracellular structure organization",
    "EMT"
  )
)

p4B <- ggplot(Clusters.Percentage, aes(x = "", y = Percentage, fill = GESs)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Paired") +
  theme_void() +
  geom_text(
    aes(label = percent(Percentage / 100)),
    size = 5, fontface = "bold",
    position = position_stack(vjust = 0.5)
  ) +
  theme(
    legend.title = element_text(size = 15, face = "bold"),
    legend.text  = element_text(size = 13, face = "bold")
  )

ggsave("Fig4B_TCGA_GESs_Pie.png", p4B, width = 8.5, height = 9, dpi = 1200)

# 6. FIGURE 4C – HALLMARK GSVA HEATMAP ----------------------------------------
cat("Generating Fig 4C Hallmark GSVA heatmap...\n")

TCGA.Seurat <- CreateSeuratObject(logCPM.all, project = "TCGA")
TCGA.Seurat <- SetIdent(TCGA.Seurat, value = clinical.w.CIBERSORTx$scRNAcluster)

TCGA.Seurat.Norm <- NormalizeData(TCGA.Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
TCGA.Average <- AverageExpression(TCGA.Seurat.Norm, return.seurat = TRUE)
TCGA.Average.table <- TCGA.Average[["RNA"]]@data

gene.sets.H <- getGeneSets(library = "H", species = "Homo sapiens")

TCGA.GSVA.Enrich <- gsva(
  as.matrix(TCGA.Average.table),
  gene.sets.H,
  min.sz = 5,
  max.sz = 500,
  kcdf  = "Gaussian"
)

write.csv(TCGA.GSVA.Enrich, "Fig4C_TCGA_GSVA_Hallmarks.csv")

ROWnames <- lapply(
  substr(gsub("_", " ", gsub("", "", rownames(TCGA.GSVA.Enrich))), 1, 500),
  function(x) bquote(bold(.(x)))
)
COLnames <- lapply(
  colnames(TCGA.GSVA.Enrich),
  function(x) bquote(bold(.(x)))
)

png("Fig4C_TCGA_Hallmarks_GESs.png", width = 10.5, height = 13, units = "in", res = 1200)
pheatmap::pheatmap(
  TCGA.GSVA.Enrich,
  fontsize     = 14,
  angle_col    = 90,
  cluster_cols = FALSE,
  fontsize_col = 18,
  labels_col   = as.expression(COLnames),
  labels_row   = as.expression(ROWnames)
)
dev.off()

# 7. FIGURE 4D – xCell CELL-TYPE ENRICHMENT HEATMAP ---------------------------
cat("Running xCell and generating Fig 4D heatmap...\n")

xCell_scores <- xCellAnalysis(logCPM.all)

xCell_filtered <- xCell_scores[!rownames(xCell_scores) %in% c(
  "MicroenvironmentScore", "StromaScore", "ImmuneScore"
), ]

xCell_filtered <- xCell_filtered[!rownames(xCell_filtered) %in% c(
  "Osteoblast", "aDC", "pDC",
  "CD4+ Tcm", "cDC", "Neurons", "CLP", "Adipocytes", "Tregs",
  "Class-switched memory B-cells", "Tgd cells", "Pericytes",
  "GMP", "pro B-cells", "iDC", "Mesangial cells",
  "Preadipocytes", "Hepatocytes", "Chondrocytes", "Th1 cells",
  "Th2 cells", "Megakaryocytes", "MPP", "Skeletal muscle", "CMP",
  "Erythrocytes", "mv Endothelial cells", "ly Endothelial cells",
  "CD4+ Tem", "Sebocytes", "Platelets", "CD8+ Tem", "CD8+ Tcm",
  "NKT"
), ]

xCell.Seurat <- CreateSeuratObject(as.matrix(xCell_filtered), project = "xCell")
xCell.Seurat <- SetIdent(xCell.Seurat, value = clinical.w.CIBERSORTx$scRNAcluster)

xCell.Average <- AverageExpression(xCell.Seurat, return.seurat = TRUE)
xCell.Average.table <- xCell.Average[["RNA"]]@data

ROWnames_x <- lapply(
  substr(gsub("_", " ", gsub("", "", rownames(xCell.Average.table))), 1, 500),
  function(x) bquote(bold(.(x)))
)
COLnames_x <- lapply(
  colnames(xCell.Average.table),
  function(x) bquote(bold(.(x)))
)

png("Fig4D_xCell_GESs_CellTypes.png", width = 7, height = 14, units = "in", res = 1200)
pheatmap::pheatmap(
  xCell.Average.table,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  angle_col    = 90,
  fontsize     = 14,
  labels_col   = as.expression(COLnames_x),
  labels_row   = as.expression(ROWnames_x)
)
dev.off()

cat("FIGURE 4 SCRIPT FINISHED (with TCGA download + processing).\n")
