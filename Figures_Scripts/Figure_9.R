# =============================================================================
# FIGURE 9: MPS_45 Independence from Clinicopathologic Features
# Panels:
#   A – Multivariate Cox: MPS_45 + AJCC M/N/T (forest plot)
#   B – Heatmap of MPS_45 gene expression across GESs
#   C – Multivariate Cox: MPS_45 + GESs (forest plot)
# =============================================================================

# 0. SETUP ---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "survival", "survminer", "forestmodel", "dplyr", "ggplot2", "Seurat", "pheatmap"
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs, ask = FALSE)
lapply(required_pkgs, library, character.only = TRUE)

options(stringsAsFactors = FALSE)
theme_set(theme_classic())
set.seed(1234)

# Assumes:
# - clinical.w.CIBERSORTx with survival + AJCC + scRNAcluster + MPS_45
# - TCGA.Voom   (genes x samples, log expression)
# - predictors  (character vector of 45 genes: MPS_45)

# 1. PREPARE CLINICAL VARIABLES AND AJCC GROUPS -------------------------------
cat("Preparing clinical variables and AJCC groups...\n")

# Survival variables (if not already defined)
if (!"overall_survival" %in% colnames(clinical.w.CIBERSORTx)) {
  clinical.w.CIBERSORTx$deceased <- clinical.w.CIBERSORTx$vital_status == "Dead"
  clinical.w.CIBERSORTx$overall_survival <- ifelse(
    clinical.w.CIBERSORTx$deceased,
    clinical.w.CIBERSORTx$days_to_death,
    clinical.w.CIBERSORTx$days_to_last_follow_up
  )
}
clinical.w.CIBERSORTx <- clinical.w.CIBERSORTx[!is.na(clinical.w.CIBERSORTx$overall_survival), ]

# AJCC M
clinical.w.CIBERSORTx$AJCC_M <- clinical.w.CIBERSORTx$ajcc_pathologic_m
clinical.w.CIBERSORTx[clinical.w.CIBERSORTx == "M1a" |
                        clinical.w.CIBERSORTx == "M1b" |
                        clinical.w.CIBERSORTx == "M1c"] <- "M1"
clinical.w.CIBERSORTx$AJCC_M[clinical.w.CIBERSORTx$AJCC_M == "M0"] <- "M0"

# AJCC N
clinical.w.CIBERSORTx$AJCC_N <- clinical.w.CIBERSORTx$ajcc_pathologic_n
clinical.w.CIBERSORTx[clinical.w.CIBERSORTx == "N1a" |
                        clinical.w.CIBERSORTx == "N1b" |
                        clinical.w.CIBERSORTx == "N1"] <- "N1"
clinical.w.CIBERSORTx[clinical.w.CIBERSORTx == "N2a" |
                        clinical.w.CIBERSORTx == "N2b" |
                        clinical.w.CIBERSORTx == "N2c" |
                        clinical.w.CIBERSORTx == "N2"] <- "N2"
clinical.w.CIBERSORTx$AJCC_N[clinical.w.CIBERSORTx$AJCC_N == "N0"] <- "N0"
clinical.w.CIBERSORTx$AJCC_N[clinical.w.CIBERSORTx$AJCC_N == "NX"] <- NA

# AJCC T
clinical.w.CIBERSORTx$AJCC_T <- clinical.w.CIBERSORTx$ajcc_pathologic_t
clinical.w.CIBERSORTx[clinical.w.CIBERSORTx == "T1a" |
                        clinical.w.CIBERSORTx == "T1b"] <- "T1"
clinical.w.CIBERSORTx[clinical.w.CIBERSORTx == "T2a" |
                        clinical.w.CIBERSORTx == "T2b" |
                        clinical.w.CIBERSORTx == "T2"] <- "T2"
clinical.w.CIBERSORTx[clinical.w.CIBERSORTx == "T3a" |
                        clinical.w.CIBERSORTx == "T3b"] <- "T3"
clinical.w.CIBERSORTx[clinical.w.CIBERSORTx == "T4a" |
                        clinical.w.CIBERSORTx == "T4b" |
                        clinical.w.CIBERSORTx == "T4"] <- "T4"
clinical.w.CIBERSORTx$AJCC_T[clinical.w.CIBERSORTx$AJCC_T %in% c("T0", "Tis", "TX")] <- NA

# Ensure factor ordering
clinical.w.CIBERSORTx$AJCC_M <- factor(clinical.w.CIBERSORTx$AJCC_M, levels = c("M0", "M1"))
clinical.w.CIBERSORTx$AJCC_N <- factor(clinical.w.CIBERSORTx$AJCC_N, levels = c("N0", "N1", "N2", "N3"))
clinical.w.CIBERSORTx$AJCC_T <- factor(clinical.w.CIBERSORTx$AJCC_T, levels = c("T1", "T2", "T3", "T4"))

# Make sure MPS_45 is factor: Low / High (as in Fig. 8)
clinical.w.CIBERSORTx$MPS_45 <- factor(clinical.w.CIBERSORTx$MPS_45, levels = c("Low", "High"))

# 2. FIGURE 9A – MULTIVARIATE COX: MPS_45 + AJCC M/N/T ------------------------
cat("Running multivariate Cox (MPS_45 + AJCC M/N/T) for Fig 9A...\n")

df_cox_A <- clinical.w.CIBERSORTx %>%
  dplyr::select(overall_survival, deceased, MPS_45, AJCC_M, AJCC_N, AJCC_T) %>%
  na.omit()

cox_A <- coxph(
  Surv(overall_survival, deceased) ~ MPS_45 + AJCC_M + AJCC_N + AJCC_T,
  data = df_cox_A
)

png("Fig9A_MPS45_Hazard_Ratio_TNM.png", width = 7.5, height = 6, units = "in", res = 1200)
forest_model(cox_A)
dev.off()

# 3. FIGURE 9B – HEATMAP OF MPS_45 GENES ACROSS GESs --------------------------

cat("Creating MPS_45 gene expression heatmap across GESs (Fig 9B)...\n")

# Build Seurat object from TCGA.Voom if not already done
if (!exists("TCGA.Seurat")) {
  TCGA.Seurat <- CreateSeuratObject(TCGA.Voom, project = "TCGA")
  TCGA.Seurat <- SetIdent(TCGA.Seurat, value = clinical.w.CIBERSORTx$scRNAcluster)
}

# Ensure same 6 GES levels
TCGA.Seurat@meta.data$scRNAcluster <- factor(
  clinical.w.CIBERSORTx$scRNAcluster,
  levels = c(
    "Anti-apoptosis",
    "Immune cell interactions",
    "Melanogenesis",
    "Ribosomal biogenesis",
    "Extracellular structure organization",
    "EMT"
  )
)
Idents(TCGA.Seurat) <- "scRNAcluster"

TCGA.Seurat <- NormalizeData(TCGA.Seurat)
TCGA.Seurat <- ScaleData(TCGA.Seurat, features = rownames(TCGA.Seurat))

# Average expression per GES
TCGA.Average <- AverageExpression(TCGA.Seurat, return.seurat = TRUE)
TCGA.Avg.mat <- TCGA.Average[["RNA"]]@data

# Keep only MPS_45 genes present
mps_genes <- predictors[predictors %in% rownames(TCGA.Avg.mat)]
TCGA.Avg.mps <- TCGA.Avg.mat[mps_genes, , drop = FALSE]

# Order genes as in your figure (optional: sort alphabetically)
TCGA.Avg.mps <- TCGA.Avg.mps[rev(rownames(TCGA.Avg.mps)), ]

png("Fig9B_MPS45_Genes_Heatmap_GESs.png", width = 6, height = 8, units = "in", res = 1200)
pheatmap::pheatmap(
  TCGA.Avg.mps,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale        = "row",
  color        = colorRampPalette(c("magenta", "black", "yellow"))(100),
  fontsize_row = 10,
  fontsize_col = 12
)
dev.off()

# 4. FIGURE 9C – MULTIVARIATE COX: MPS_45 + GESs ------------------------------
cat("Running multivariate Cox (MPS_45 + GESs) for Fig 9C...\n")

# Define GES factor (same ordering as other figures)
clinical.w.CIBERSORTx$GES <- factor(
  clinical.w.CIBERSORTx$scRNAcluster,
  levels = c(
    "Immune cell interactions",
    "Anti-apoptosis",
    "Melanogenesis",
    "Ribosomal biogenesis",
    "Extracellular structure organization",
    "EMT"
  ),
  ordered = TRUE
)

df_cox_C <- clinical.w.CIBERSORTx %>%
  dplyr::select(overall_survival, deceased, MPS_45, GES) %>%
  na.omit()

cox_C <- coxph(
  Surv(overall_survival, deceased) ~ MPS_45 + GES,
  data = df_cox_C
)

png("Fig9C_MPS45_Hazard_Ratio_with_GESs.png", width = 7.5, height = 6, units = "in", res = 1200)
forest_model(cox_C)
dev.off()

cat("FIGURE 9 SCRIPT FINISHED.\n")
cat("Generated:\n")
cat("- Fig9A_MPS45_Hazard_Ratio_TNM.png\n")
cat("- Fig9B_MPS45_Genes_Heatmap_GESs.png\n")
cat("- Fig9C_MPS45_Hazard_Ratio_with_GESs.png\n")
