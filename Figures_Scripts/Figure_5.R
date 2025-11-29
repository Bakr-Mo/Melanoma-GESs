# =============================================================================
# FIGURE 5: Overlap of Melanoma GESs with Published Signatures
# Panels:
#   A – Heatmaps: immune and high-immune signatures across TCGA GESs
#   B – Heatmaps: MITF-low and proliferative signatures across TCGA GESs
#   C – MITF dot plot across GESs
#   D – UpSet plot of patient-level overlap between our GESs and TCGA 2015 clusters
# =============================================================================

# 0. SETUP ---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "Seurat", "dplyr", "ggplot2", "data.table",
  "UpSetR", "ComplexUpset"
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs, ask = FALSE)
lapply(required_pkgs, library, character.only = TRUE)

options(stringsAsFactors = FALSE)
theme_set(theme_classic())
set.seed(1234)

# Assumes:
# - TCGA.Seurat and TCGA.Average are in memory from Figure 4
#   (TCGA.Average <- AverageExpression(TCGA.Seurat.Norm, return.seurat = TRUE))
# - clinical.w.CIBERSORTx exists with $scRNAcluster labels
# - MITF.GSE115978 exists (Seurat object with scRNA melanoma clusters)

# 1. LOAD PUBLISHED SIGNATURES (TCGA 2015, JONSSON 2010) ----------------------
cat("Loading external GES signature tables...\n")

TCGA2015 <- read.csv("TCGA2015.csv", header = TRUE, sep = ",")    # contains gene lists (e.g. MITF.low, Immune, etc.)
Jonsson2010 <- read.csv("Jonsson2010.csv", header = TRUE, sep = ",")

# 2. PREPARE TCGA AVERAGE EXPRESSION OBJECT -----------------------------------
cat("Preparing TCGA average expression object...\n")

TCGA.Seurat.genes <- rownames(TCGA.Seurat)
TCGA.Seurat <- ScaleData(TCGA.Seurat, features = TCGA.Seurat.genes)

# Ensure GES order matches manuscript
levels(TCGA.Average) <- c(
  "Anti-apoptosis",
  "Immune cell interactions",
  "Melanogenesis",
  "Ribosomal biogenesis",
  "Extracellular structure organization",
  "EMT"
)

# 3. FIGURE 5A/5B – HEATMAPS OF PUBLISHED SIGNATURES -------------------------
# You need to define 'predictors' as union of representative genes, or use
# columns like TCGA2015$Immune, TCGA2015$HighImmune, TCGA2015$MITF.low, etc.

# Example: 'predictors' = all signature genes from TCGA2015 and Jonsson2010
# (edit to exactly match the two panels A and B)
predictors <- unique(na.omit(c(
  TCGA2015$immune,       # replace with actual column names
  TCGA2015$high_immune,  # e.g. 'HighImmune' if present
  TCGA2015$MITF.low,
  Jonsson2010$proliferative
)))

predictors <- predictors[predictors %in% rownames(TCGA.Average[["RNA"]]@data)]

cat("Plotting signature heatmaps (Fig 5A/5B)...\n")

p5AB <- DoHeatmap(
  TCGA.Average,
  features   = predictors,
  angle      = 90,
  label      = FALSE,
  draw.lines = FALSE
) +
  theme(
    axis.text.y = element_text(size = 8, colour = "black", face = "bold.italic"),
    axis.text.x = element_text(size = 14, colour = "black", face = "bold")
  )

ggsave("Fig5AB_SignatureHeatmaps_TCGA_GESs.png", p5AB, width = 7, height = 12, dpi = 1200)

# If you want separate panels for immune/high-immune vs MITF-low/proliferative:
# immune_predictors     <- na.omit(TCGA2015$immune)
# highimmune_predictors <- na.omit(TCGA2015$high_immune)
# MITFlow_predictors    <- na.omit(TCGA2015$MITF.low)
# prolif_predictors     <- na.omit(Jonsson2010$proliferative)
# and then run DoHeatmap() twice with those sets.

# 4. FIGURE 5C – MITF DOT PLOT ACROSS TCGA GESs -------------------------------
cat("Plotting MITF dot plot across GESs (Fig 5C)...\n")

p5C <- DotPlot(
  TCGA.Average,
  features   = "MITF",
  dot.scale  = 10,
  cols       = c("blue", "red")
) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(size = 25, colour = "black", face = "bold.italic"),
    axis.text.y  = element_text(size = 20, colour = "black", face = "bold"),
    legend.text  = element_text(size = 15),
    legend.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1))

ggsave("Fig5C_MITF_DotPlot_TCGA_GESs.png", p5C, width = 5, height = 4, dpi = 1200)

# Optional AXL in scRNA melanoma clusters (not directly used in Fig 5 but useful)
p5C_scRNA <- DotPlot(
  MITF.GSE115978,
  features   = "AXL",
  dot.scale  = 10,
  cols       = c("blue", "red")
) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(size = 25, colour = "black", face = "bold.italic"),
    axis.text.y  = element_text(size = 20, colour = "black", face = "bold"),
    legend.text  = element_text(size = 15),
    legend.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1))

ggsave("Fig5C_AXL_DotPlot_scRNA_GESs.png", p5C_scRNA, width = 5, height = 4, dpi = 1200)

# 5. FIGURE 5A/B – HEATMAPS USING TCGA 2015 PATIENT SUBSET --------------------
cat("Preparing TCGA 2015 patient subset for heatmaps...\n")

TCGA2015.Patients <- read.csv("TCGA2015_Patients.csv", header = TRUE, sep = ",")
TCGA2015.Patients.329 <- TCGA2015.Patients[
  TCGA2015.Patients$RNASEQ.CLUSTER_CONSENHIER != "-", ]

TCGA.Clinical <- clinical.w.CIBERSORTx
TCGA.Clinical$ID <- gsub(".{1}$", "", TCGA.Clinical$sample)

ID.329 <- TCGA.Clinical[TCGA.Clinical$ID %in% TCGA2015.Patients.329$Name, ]
TCGA.Voom <- logCPM.all   # from Figure 4 pipeline

TCGA.Voom.329 <- TCGA.Voom[, colnames(TCGA.Voom) %in% rownames(ID.329)]

TCGA.Voom.329.Seurat <- CreateSeuratObject(TCGA.Voom.329, project = "TCGA_329")
TCGA.Voom.329.Seurat <- SetIdent(TCGA.Voom.329.Seurat, value = ID.329$scRNAcluster)

TCGA.Voom.329.Seurat.Average <- AverageExpression(TCGA.Voom.329.Seurat, return.seurat = TRUE)

levels(TCGA.Voom.329.Seurat.Average) <- c(
  "Anti-apoptosis",
  "Immune cell interactions",
  "Melanogenesis",
  "Ribosomal biogenesis",
  "Extracellular structure organization",
  "EMT"
)

# Example: show MITF-low signature across our GESs (panel B right heatmap)
MITFlow_genes <- na.omit(TCGA2015$MITF.low)
MITFlow_genes <- MITFlow_genes[MITFlow_genes %in% rownames(TCGA.Voom.329.Seurat.Average[["RNA"]]@data)]

p5B_MITF <- DoHeatmap(
  TCGA.Voom.329.Seurat.Average,
  features   = MITFlow_genes,
  angle      = 90,
  label      = FALSE,
  draw.lines = FALSE
) +
  theme(axis.text.y = element_blank())

ggsave("Fig5B_MITFlow_TCGA2015_vs_GESs.png", p5B_MITF, width = 4, height = 8, dpi = 1200)

# Similar blocks can be added for “immune” / “high immune” signatures.

# 6. FIGURE 5D – UPSET PLOT OF PATIENT OVERLAP --------------------------------
cat("Generating UpSet plot for patient overlap (Fig 5D)...\n")

Keratin      <- TCGA2015.Patients.329$Name[
  TCGA2015.Patients.329$RNASEQ.CLUSTER_CONSENHIER == "keratin"
]
Immune2015   <- TCGA2015.Patients.329$Name[
  TCGA2015.Patients.329$RNASEQ.CLUSTER_CONSENHIER == "immune"
]
MITF.low2015 <- TCGA2015.Patients.329$Name[
  TCGA2015.Patients.329$RNASEQ.CLUSTER_CONSENHIER == "MITF-low"
]

Anti.apoptosis <- ID.329$ID[ID.329$scRNAcluster == "Anti-apoptosis"]
Immune.cells.interactions <- ID.329$ID[ID.329$scRNAcluster == "Immune cell interactions"]
Melanogenesis <- ID.329$ID[ID.329$scRNAcluster == "Melanogenesis"]
Ribosomal.biogenesis <- ID.329$ID[ID.329$scRNAcluster == "Ribosomal biogenesis"]
Extracellular.structure.organization <- ID.329$ID[ID.329$scRNAcluster == "Extracellular structure organization"]
EMT <- ID.329$ID[ID.329$scRNAcluster == "EMT"]

x_list <- list(
  Keratin      = Keratin,
  Immune       = Immune2015,
  MITF.low     = MITF.low2015,
  Anti.apoptosis = Anti.apoptosis,
  `Immune cell interactions` = Immune.cells.interactions,
  Melanogenesis = Melanogenesis,
  `Ribosomal biogenesis` = Ribosomal.biogenesis,
  `Extracellular structure organization` = Extracellular.structure.organization,
  EMT = EMT
)

names(x_list) <- c(
  "keratin (TCGA, 2015)",
  "immune (TCGA, 2015)",
  "MITF-low (TCGA, 2015)",
  "Anti-apoptosis",
  "Immune cell interactions",
  "Melanogenesis",
  "Ribosomal biogenesis",
  "Extracellular structure organization",
  "EMT"
)

upset_list <- fromList(x_list)

png("Fig5D_UpSet_TCGA2015_vs_GESs.png", width = 9, height = 6, units = "in", res = 1200)
ComplexUpset::upset(
  upset_list,
  c(
    "EMT",
    "Extracellular structure organization",
    "Ribosomal biogenesis",
    "Melanogenesis",
    "Immune cell interactions",
    "Anti-apoptosis",
    "MITF-low (TCGA, 2015)",
    "keratin (TCGA, 2015)",
    "immune (TCGA, 2015)"
  ),
  min_size          = 5,
  sort_intersections = FALSE,
  width_ratio       = 0.2,
  base_annotations  = list(
    "Intersection ratio" = intersection_ratio(
      text_mapping = aes(label = !!upset_text_percentage())
    )
  ),
  sort_sets         = FALSE
)
dev.off()

cat("FIGURE 5 SCRIPT FINISHED.\n")
cat("Generated panels:\n")
cat("- Fig5AB_SignatureHeatmaps_TCGA_GESs.png (and optional per-signature heatmaps)\n")
cat("- Fig5C_MITF_DotPlot_TCGA_GESs.png (plus AXL DotPlot if desired)\n")
cat("- Fig5B_MITFlow_TCGA2015_vs_GESs.png (TCGA 2015 MITF-low vs our GESs)\n")
cat("- Fig5D_UpSet_TCGA2015_vs_GESs.png\n")
