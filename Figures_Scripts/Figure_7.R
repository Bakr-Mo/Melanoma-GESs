# =============================================================================
# FIGURE 7: Univariate Cox Regression for Ribosomal Biogenesis DEGs
# Goal:
#   - Fit univariate Cox models for each Ribo-Only scRNA DEG in TCGA.
#   - Select genes with |log2FC| > 1 and p < 0.05.
#   - Plot HRs and 95% CIs as a forest plot (colored by p-value).
# =============================================================================

# 0. SETUP ---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "survival", "survminer", "plyr", "edgeR", "dplyr", "ggplot2"
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs, ask = FALSE)
lapply(required_pkgs, library, character.only = TRUE)

options(stringsAsFactors = FALSE)
theme_set(theme_classic())
set.seed(1234)

# Assumes:
# - all.genes.pc (Ensembl + gene names; protein coding)
# - all.gene.names (external gene symbols matching all.genes.pc rows)
# - clinical (from TCGA SKCM SummarizedExperiment)
# - clinical.w.CIBERSORTx (same clinical plus GES labels)
# - Ribo.Only.scRNA.csv with column 'Markers' listing Ribo DEGs

# 1. LOAD RIBOSOMAL-BIOGENESIS MARKERS ----------------------------------------
cat("Loading ribosomal-biogenesis marker lists...\n")

Ribo.Only.scRNA <- read.csv("Ribo_Markers_scRNA_logfc1.csv", header = TRUE, sep = ",")
ribo_markers <- unique(na.omit(Ribo.Only.scRNA$Markers))

# 2. BUILD TCGA.Voom EXPRESSION (EDGE R + VOOM) -------------------------------
cat("Building TCGA.Voom expression matrix (TMM + voom)...\n")

Not.Norm.TCGA <- all.genes.pc
rownames(Not.Norm.TCGA) <- all.gene.names

sample.type <- clinical$sample_type
DGE <- DGEList(Not.Norm.TCGA, group = sample.type)
group <- factor(sample.type)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

keep <- filterByExpr(DGE, design)
DGE <- DGE[keep, , keep.lib.sizes = FALSE]

DGE.normalized <- calcNormFactors(DGE, method = "TMM")
Voom <- voom(DGE.normalized, design)
TCGA.Voom <- Voom$E  # genes x samples (log2 CPM-ish)

# 3. PREPARE BINARY UP/DOWN MATRIX FOR RIBO MARKERS ---------------------------
cat("Transforming TCGA.Voom into up/down calls for ribo markers...\n")

Trans.TCGA.Voom <- t(TCGA.Voom)  # samples x genes

# Median-based up/down per gene
Trans.TCGA.Voom.upANDdown <- apply(
  Trans.TCGA.Voom,
  2,
  function(x) ifelse(x > median(x), "up", "down")
) %>%
  as.data.frame()

# Restrict to ribosomal biogenesis markers (intersection with TCGA genes)
ribo_markers_in_tcga <- colnames(Trans.TCGA.Voom.upANDdown)[
  colnames(Trans.TCGA.Voom.upANDdown) %in% ribo_markers
]

Trans.TCGA.Voom.upANDdown.Ribo.Only.scRNA <- Trans.TCGA.Voom.upANDdown[, ribo_markers_in_tcga, drop = FALSE]

# Make sure column names are syntactically valid (for HLA and others)
colnames(Trans.TCGA.Voom.upANDdown.Ribo.Only.scRNA) <- make.names(
  colnames(Trans.TCGA.Voom.upANDdown.Ribo.Only.scRNA)
)

# 4. MERGE EXPRESSION WITH CLINICAL + SURVIVAL --------------------------------
cat("Merging expression with clinical and defining survival...\n")

clinical.w.CIBERSORTx$deceased <- clinical.w.CIBERSORTx$vital_status == "Dead"
clinical.w.CIBERSORTx$overall_survival <- ifelse(
  clinical.w.CIBERSORTx$deceased,
  clinical.w.CIBERSORTx$days_to_death,
  clinical.w.CIBERSORTx$days_to_last_follow_up
)

clinical.w.CIBERSORTx <- clinical.w.CIBERSORTx[!is.na(clinical.w.CIBERSORTx$overall_survival), ]

# Align samples between clinical and expression
common_samples <- intersect(
  rownames(clinical.w.CIBERSORTx),
  rownames(Trans.TCGA.Voom.upANDdown.Ribo.Only.scRNA)
)

clinical_cox <- clinical.w.CIBERSORTx[common_samples, ]
expr_cox     <- Trans.TCGA.Voom.upANDdown.Ribo.Only.scRNA[common_samples, , drop = FALSE]

clinical.w.CIBERSORTx.Ribo.Only.scRNA <- cbind(clinical_cox, expr_cox)

# 5. UNIVARIATE COX REGRESSION PER GENE ---------------------------------------
cat("Running univariate Cox regression for ribo markers...\n")

gene_cols <- colnames(expr_cox)

formulas <- sapply(
  gene_cols,
  USE.NAMES = TRUE,
  function(x) as.formula(paste("Surv(overall_survival, deceased) ~", x))
)

models <- lapply(
  formulas,
  function(fm) coxph(fm, data = clinical.w.CIBERSORTx.Ribo.Only.scRNA)
)

res <- lapply(
  models,
  function(m) {
    cbind(
      HR    = exp(coef(m)),
      exp(confint(m)),
      Pval  = coef(summary(m))[5]
    )
  }
)

HR.Table.Ribo.Only.scRNA <- plyr::ldply(res, data.frame)
colnames(HR.Table.Ribo.Only.scRNA) <- c("Gene", "HR", "lower", "upper", "Pval")

# Filter to significant genes (p < 0.05)
HR.Table.Ribo.Only.scRNA <- HR.Table.Ribo.Only.scRNA[HR.Table.Ribo.Only.scRNA$Pval < 0.05, ]

# Optional: order by HR or Pval
HR.Table.Ribo.Only.scRNA <- HR.Table.Ribo.Only.scRNA[order(HR.Table.Ribo.Only.scRNA$HR), ]

write.csv(HR.Table.Ribo.Only.scRNA, "Fig7_HR_Table_Ribo_Only_scRNA.csv", row.names = FALSE)

# 6. FOREST PLOT (FIGURE 7) ---------------------------------------------------
cat("Plotting forest plot (Fig 7)...\n")

p7 <- ggplot(HR.Table.Ribo.Only.scRNA, aes(x = HR, y = Gene)) +
  geom_pointrange(
    aes(xmin = lower, xmax = upper, colour = Pval),
    position = position_dodge(width = 0.5),
    size     = 2,
    fatten   = 1.5
  ) +
  geom_vline(xintercept = 1, linetype = 5) +
  scale_color_gradient(
    name  = "Pval",
    low   = "#084081",
    high  = "#7bccc4",
    limits = c(0.01, 0.04) # adjust to emphasise 0.01–0.04 as in the legend
  ) +
  theme_classic() +
  theme(
    strip.placement   = "outside",
    strip.text        = element_text(angle = 180),
    strip.background  = element_blank(),
    panel.spacing     = unit(1, "mm"),
    axis.title.x      = element_text(size = 25, colour = "black", face = "bold"),
    axis.text.x       = element_text(size = 20, colour = "black", face = "bold"),
    axis.text.y       = element_text(size = 15, colour = "black", face = "bold.italic"),
    legend.text       = element_text(size = 13, colour = "black", face = "bold"),
    legend.title      = element_text(size = 20, face = "bold")
  ) +
  xlab("HR") +
  ylab("")

ggsave("Fig7_HR_Ribo_LogFC1_scRNA.png", p7, width = 8, height = 11, dpi = 1200)

cat("FIGURE 7 SCRIPT FINISHED.\n")
cat("Generated:\n")
cat("- Fig7_HR_Ribo_LogFC1_scRNA.png\n")
cat("- Fig7_HR_Table_Ribo_Only_scRNA.csv\n")
