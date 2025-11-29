# =============================================================================
# FIGURE 8: Prognostic Signature MPS_45 in TCGA-SKCM and External Cohorts
# Panels:
#   A – Ridge (α=0) Cox regression for 45 ribosomal-bio DEGs (MPS_45)
#   B – KM: TCGA-SKCM High vs Low MPS_45
#   C – KM: GSE65904 High vs Low MPS_45
#   D – KM: GSE19234 High vs Low MPS_45
#   E – KM: GSE53118 High vs Low MPS_45
#   F – Time-dependent ROC for TCGA (1,3,5 years)
# =============================================================================

# 0. SETUP ---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "glmnet", "survival", "survminer", "edgeR", "dplyr", "ggplot2",
  "pROC", "singscore", "GEOquery", "Biobase", "forcats"
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs, ask = FALSE)
lapply(required_pkgs, library, character.only = TRUE)

options(stringsAsFactors = FALSE)
theme_set(theme_classic())
set.seed(1234)

# Assumes:
# - TCGA.Voom      : genes x samples, log2 expression from voom
# - clinical.w.CIBERSORTx: includes vital_status, days_to_death, days_to_last_follow_up
# - HR.Table.Ribo.Only.scRNA: from Figure 7, with columns Gene, HR, Pval
# - predictors: character vector of 45 genes selected from ridge (MPS_45)

# 1. RISK MODEL CONSTRUCTION IN TCGA (RIDGE COX, FIG 8A) ----------------------

cat("Constructing ridge Cox model to define MPS_45...\n")

# Use prognostic ribo genes with HR > 1.3 and p < 0.05
HR.Table.Ribo.Only.scRNA.HR <- HR.Table.Ribo.Only.scRNA[
  HR.Table.Ribo.Only.scRNA$Pval < 0.05 & HR.Table.Ribo.Only.scRNA$HR > 1.3, ]
Ribo.Genes <- HR.Table.Ribo.Only.scRNA.HR$Gene

# Extract expression for these genes from TCGA.Voom
Ribo.Genes.Expression <- TCGA.Voom[rownames(TCGA.Voom) %in% Ribo.Genes, , drop = FALSE]
Ribo.Genes.Expression <- as.data.frame(t(Ribo.Genes.Expression))  # samples x genes

# Survival table
TCGA_clin <- data.frame(
  time   = clinical.w.CIBERSORTx$overall_survival,
  status = ifelse(clinical.w.CIBERSORTx$vital_status == "Dead", 1, 0)
)
rownames(TCGA_clin) <- rownames(clinical.w.CIBERSORTx)
TCGA_clin$time[TCGA_clin$time <= 0] <- NA
TCGA_clin <- na.omit(TCGA_clin)

# Align expression with survival table
Ribo.Genes.Expression.Adj <- Ribo.Genes.Expression[
  rownames(Ribo.Genes.Expression) %in% rownames(TCGA_clin),
]

# Ridge Cox
x <- as.matrix(Ribo.Genes.Expression.Adj)
y <- Surv(TCGA_clin$time[rownames(TCGA_clin) %in% rownames(Ribo.Genes.Expression.Adj)],
          TCGA_clin$status[rownames(TCGA_clin) %in% rownames(Ribo.Genes.Expression.Adj)])

fit <- glmnet(x, y, family = "cox", alpha = 0)

# Plot lambda path (not saved)
# plot(fit)

cv.fit <- cv.glmnet(x, y, family = "cox", alpha = 0)

png("Fig8A_Ridge_Regression_Ribo_scRNA.png", width = 7, height = 6, units = "in", res = 1200)
plot(cv.fit, cex.lab = 1.5, cex.axis = 1.2, font.lab = 2)
dev.off()

idealLambda <- cv.fit$lambda.min
co <- coef(cv.fit, s = idealLambda)
predictors <- rownames(co)[co != 0]  # MPS_45 genes
write.table(predictors, "MPS45_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 2. DEFINE MPS_45 SCORE AND GROUPS IN TCGA (FIG 8B, 8F) ----------------------

cat("Scoring MPS_45 in TCGA...\n")

rankData <- rankGenes(TCGA.Voom)
score_tcga <- simpleScore(rankData, upSet = predictors)

score_groups <- apply(score_tcga, 2, function(x) ifelse(x > median(x), "High", "Low")) %>%
  as.data.frame()

clinical.w.CIBERSORTx$MPS_45 <- ifelse(score_groups$TotalScore == "High", "High", "Low")

# Prepare survival variables (if not already done)
clinical.w.CIBERSORTx$deceased <- clinical.w.CIBERSORTx$vital_status == "Dead"
clinical.w.CIBERSORTx$overall_survival <- ifelse(
  clinical.w.CIBERSORTx$deceased,
  clinical.w.CIBERSORTx$days_to_death,
  clinical.w.CIBERSORTx$days_to_last_follow_up
)

clinical.w.CIBERSORTx <- clinical.w.CIBERSORTx[!is.na(clinical.w.CIBERSORTx$overall_survival), ]

# Reverse for Cox so High is reference, as in your code
clinical.w.CIBERSORTx$MPS_45  <- clinical.w.CIBERSORTx$MPS_45 %>%
  as.factor() %>% forcats::fct_rev()

# KM plot (Fig 8B)
fit_tcga_mps <- survfit(Surv(overall_survival, deceased) ~ MPS_45,
                        data = clinical.w.CIBERSORTx)
pval_tcga   <- surv_pvalue(fit_tcga_mps, data = clinical.w.CIBERSORTx)$pval
print(pval_tcga)

p_tcga <- ggsurvplot(
  fit_tcga_mps,
  data       = clinical.w.CIBERSORTx,
  pval       = FALSE,
  risk.table = TRUE,
  size       = 2.5,
  legend.labs  = c("Low MPS_45", "High MPS_45"),
  legend.title = "",
  font.legend  = c(25, "bold"),
  font.x       = c(25,  "bold"),
  font.y       = c(25,  "bold"),
  legend       = c(0.7, 0.85),
  palette      = c("blue", "red")
)

ggsave(
  filename = "Fig8B_MPS45_KM_TCGA.png",
  plot     = print(p_tcga$plot, newpage = FALSE),
  device   = "png",
  width    = 6,
  height   = 5,
  dpi      = 1200
)

# 3. TIME-DEPENDENT ROC IN TCGA (FIG 8F) --------------------------------------
cat("Computing time-dependent ROC (approximate via logistic models at 1,3,5y)...\n")

# Define 1,3,5-year alive/dead labels
AliveDead <- clinical.w.CIBERSORTx[, c("days_to_death", "days_to_last_follow_up", "vital_status")]
cut_years <- c("1y" = 365, "3y" = 3 * 365, "5y" = 5 * 365)

make_status_at <- function(days_cut) {
  dead <- AliveDead$vital_status == "Dead" & AliveDead$days_to_death < days_cut
  alive <- AliveDead$vital_status == "Alive" & AliveDead$days_to_last_follow_up > days_cut
  factor(ifelse(dead, 1, ifelse(alive, 0, NA)))
}

status_1y <- make_status_at(cut_years["1y"])
status_3y <- make_status_at(cut_years["3y"])
status_5y <- make_status_at(cut_years["5y"])

# Logistic models with MPS_45 score (TotalScore) as predictor
mps_vector <- score_tcga$TotalScore[rownames(clinical.w.CIBERSORTx)]

df_1y <- data.frame(status = as.numeric(as.character(status_1y)), MPS = mps_vector)
df_1y <- df_1y[!is.na(df_1y$status), ]
glm1 <- glm(status ~ MPS, data = df_1y, family = binomial)
roc1 <- roc(glm1$y, fitted(glm1), smooth = FALSE)

df_3y <- data.frame(status = as.numeric(as.character(status_3y)), MPS = mps_vector)
df_3y <- df_3y[!is.na(df_3y$status), ]
glm3 <- glm(status ~ MPS, data = df_3y, family = binomial)
roc3 <- roc(glm3$y, fitted(glm3), smooth = FALSE)

df_5y <- data.frame(status = as.numeric(as.character(status_5y)), MPS = mps_vector)
df_5y <- df_5y[!is.na(df_5y$status), ]
glm5 <- glm(status ~ MPS, data = df_5y, family = binomial)
roc5 <- roc(glm5$y, fitted(glm5), smooth = FALSE)

png("Fig8F_ROC_TCGA_MPS45.png", width = 6, height = 6, units = "in", res = 1200)
par(pty = "s")
plot.roc(
  roc1,
  legacy.axes = TRUE,
  xlab = "False Positive Rate",
  ylab = "True Positive Rate",
  col  = "red",
  lwd  = 5,
  print.auc     = TRUE,
  cex.lab       = 2,
  cex.axis      = 1.1,
  font.lab      = 2,
  font.axis     = 2,
  print.auc.cex = 2
)
plot.roc(
  roc3,
  col  = "blue",
  lwd  = 5,
  add  = TRUE,
  print.auc     = TRUE,
  print.auc.adj = c(0, 4),
  print.auc.cex = 2
)
plot.roc(
  roc5,
  col  = "darkgreen",
  lwd  = 5,
  add  = TRUE,
  print.auc     = TRUE,
  print.auc.adj = c(0, 2.5),
  print.auc.cex = 2
)

legend(
  "bottomright",
  legend = c("1 year", "3 years", "5 years"),
  col    = c("red", "blue", "darkgreen"),
  lwd    = 6,
  cex    = 1.2
)
dev.off()

# 4. VALIDATION IN GSE65904 (FIG 8C) ------------------------------------------
cat("Validating MPS_45 in GSE65904...\n")

GSE65904 <- getGEO("GSE65904", GSEMatrix = TRUE)
gse <- GSE65904[[1]]

expr <- Biobase::exprs(gse)
matrix_65904 <- edgeR::cpm(expr, prior.count = 2, log = TRUE)

annot_65904 <- Biobase::fData(gse)
annotation_65904 <- dplyr::select(annot_65904, "Symbol")
rownames(matrix_65904) <- annotation_65904$Symbol

clin_65904 <- Biobase::pData(gse)
clin_65904_surv <- clin_65904[, c(
  "disease specific survival (1=death, 0=alive):ch1",
  "disease specific survival in days:ch1"
)]
colnames(clin_65904_surv) <- c("deceased", "overall_survival")
clin_65904_surv <- clin_65904_surv %>%
  mutate(
    deceased         = as.numeric(deceased),
    overall_survival = as.numeric(overall_survival)
  )

rankData_65904 <- rankGenes(matrix_65904)
score_65904 <- simpleScore(rankData_65904, upSet = predictors)
score_groups_65904 <- apply(score_65904, 2, function(x) ifelse(x > median(x), "High", "Low")) %>%
  as.data.frame()

clin_65904_surv$MPS_45 <- ifelse(score_groups_65904$TotalScore == "High", "High", "Low")

fit_65904 <- survfit(Surv(overall_survival, deceased) ~ MPS_45,
                     data = clin_65904_surv)
pval_65904 <- surv_pvalue(fit_65904, data = clin_65904_surv)$pval
print(pval_65904)

p_65904 <- ggsurvplot(
  fit_65904,
  data       = clin_65904_surv,
  pval       = FALSE,
  risk.table = TRUE,
  size       = 2.5,
  legend.labs  = c("High MPS_45", "Low MPS_45"),
  legend.title = "",
  font.legend  = c(25, "bold"),
  font.x       = c(25, "bold"),
  font.y       = c(25, "bold"),
  legend       = c(0.65, 0.9),
  palette      = c("red", "blue")
)

ggsave(
  filename = "Fig8C_MPS45_KM_GSE65904.png",
  plot     = print(p_65904$plot, newpage = FALSE),
  device   = "png",
  width    = 6,
  height   = 5,
  dpi      = 1200
)

# 5. VALIDATION IN GSE19234 (FIG 8D) ------------------------------------------
cat("Validating MPS_45 in GSE19234...\n")

GSE19234 <- getGEO("GSE19234", GSEMatrix = TRUE)
gse <- GSE19234[[1]]

expr <- Biobase::exprs(gse)
matrix_19234 <- edgeR::cpm(expr, prior.count = 2, log = TRUE)

annot_19234 <- Biobase::fData(gse)
annotation_19234 <- dplyr::select(annot_19234, "Gene Symbol")
rownames(matrix_19234) <- annotation_19234$`Gene Symbol`

clin_19234 <- Biobase::pData(gse)
clin_19234_surv <- clin_19234[, c(
  "staus dead or alive:ch1",
  "days lived since metastasis:ch1"
)]
colnames(clin_19234_surv) <- c("deceased", "overall_survival")
clin_19234_surv <- clin_19234_surv %>%
  mutate(
    deceased         = as.numeric(deceased),
    overall_survival = as.numeric(overall_survival)
  )

rankData_19234 <- rankGenes(matrix_19234)
score_19234 <- simpleScore(rankData_19234, upSet = predictors)
score_groups_19234 <- apply(score_19234, 2, function(x) ifelse(x > median(x), "High", "Low")) %>%
  as.data.frame()

clin_19234_surv$MPS_45 <- ifelse(score_groups_19234$TotalScore == "High", "High", "Low")

fit_19234 <- survfit(Surv(overall_survival, deceased) ~ MPS_45,
                     data = clin_19234_surv)
pval_19234 <- surv_pvalue(fit_19234, data = clin_19234_surv)$pval
print(pval_19234)

p_19234 <- ggsurvplot(
  fit_19234,
  data       = clin_19234_surv,
  pval       = FALSE,
  risk.table = TRUE,
  size       = 2.5,
  legend.labs  = c("Low MPS_45", "High MPS_45"),
  legend.title = "",
  font.legend  = c(25, "bold"),
  font.x       = c(25, "bold"),
  font.y       = c(25, "bold"),
  legend       = c(0.65, 0.9),
  palette      = c("blue", "red")
)

ggsave(
  filename = "Fig8D_MPS45_KM_GSE19234.png",
  plot     = print(p_19234$plot, newpage = FALSE),
  device   = "png",
  width    = 6,
  height   = 5,
  dpi      = 1200
)

# 6. VALIDATION IN GSE53118 (FIG 8E) ------------------------------------------
cat("Validating MPS_45 in GSE53118...\n")

GSE53118 <- getGEO("GSE53118", GSEMatrix = TRUE)
gse <- GSE53118[[1]]

expr <- Biobase::exprs(gse)
matrix_53118 <- edgeR::cpm(expr, prior.count = 2, log = TRUE)

annot_53118 <-
  