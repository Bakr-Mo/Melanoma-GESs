# =============================================================================
# FIGURE 6: Prognostic Impact of Melanoma GESs (TCGA-SKCM)
# Panels:
#   A – KM curves for all 6 GESs (global log-rank)
#   B – Bar chart of Alive >4y vs Dead <4y proportions by GES
#   C–H – KM curves for each GES vs rest of cohort
# =============================================================================

# 0. SETUP ---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "survival", "survminer", "dplyr", "ggplot2", "tidyr"
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs, ask = FALSE)
lapply(required_pkgs, library, character.only = TRUE)

options(stringsAsFactors = FALSE)
theme_set(theme_classic())
set.seed(1234)

# Assumes clinical.w.CIBERSORTx is already in memory from the TCGA script.
# If not, load it from an .rds or .csv here:
# clinical.w.CIBERSORTx <- readRDS("clinical_w_CIBERSORTx_TCGA.rds")

# 1. DEFINE SURVIVAL VARIABLES -------------------------------------------------
cat("Defining survival variables...\n")

clinical.w.CIBERSORTx$deceased <- clinical.w.CIBERSORTx$vital_status == "Dead"

clinical.w.CIBERSORTx$overall_survival <- ifelse(
  clinical.w.CIBERSORTx$deceased,
  clinical.w.CIBERSORTx$days_to_death,
  clinical.w.CIBERSORTx$days_to_last_follow_up
)

# Drop samples without survival time
clinical.w.CIBERSORTx <- clinical.w.CIBERSORTx[!is.na(clinical.w.CIBERSORTx$overall_survival), ]

# Ensure GES factor order
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

# 2. FIGURE 6A – KM CURVES FOR ALL 6 GESs -------------------------------------
cat("Fitting global GES survival model (Fig 6A)...\n")

fit_all <- survfit(
  Surv(overall_survival, deceased) ~ scRNAcluster,
  data = clinical.w.CIBERSORTx
)

pval_all <- surv_pvalue(fit_all, data = clinical.w.CIBERSORTx)$pval
print(pval_all)

surv_plot_all <- ggsurvplot(
  fit_all,
  data       = clinical.w.CIBERSORTx,
  pval       = FALSE,
  risk.table = TRUE,
  size       = 3,
  legend.labs = c(
    "Anti-apoptosis",
    "Immune cell interactions",
    "Melanogenesis",
    "Ribosomal biogenesis",
    "Extracellular structure organization",
    "EMT"
  ),
  legend.title  = "GESs",
  font.legend   = c(25, "bold"),
  font.x        = c(30, "bold"),
  font.y        = c(30, "bold"),
  legend        = c(0.67, 0.83),
  palette       = "lancet"
)

surv_plot_all$plot <- surv_plot_all$plot +
  annotate(
    "text",
    x      = 2600,
    y      = 0.1,
    label  = "p < 0.0001",
    cex    = 10,
    vjust  = 0,
    hjust  = 1.1,
    fontface = 4,
    color  = "blue"
  ) +
  xlab("Time") +
  ylab("Survival probability")

ggsave(
  filename = "Fig6A_GESs_Global_Survival.png",
  plot     = print(surv_plot_all$plot, newpage = FALSE),
  device   = "png",
  width    = 11,
  height   = 8,
  dpi      = 1200
)

# 3. ALIVE/DEAD WITHIN 4 YEARS (FOR PANEL B) ----------------------------------
cat("Deriving 4-year alive/dead classifications...\n")

AliveAndDead <- clinical.w.CIBERSORTx[, c("days_to_death", "days_to_last_follow_up", "vital_status", "scRNAcluster")]

AliveAndDead$deadStatus <- ifelse(
  AliveAndDead$vital_status == "Dead" &
    !is.na(AliveAndDead$days_to_death) &
    AliveAndDead$days_to_death < 1460,
  "Dead", "Not yet"
)

AliveAndDead$AliveStatus <- ifelse(
  AliveAndDead$vital_status == "Alive" &
    !is.na(AliveAndDead$days_to_last_follow_up) &
    AliveAndDead$days_to_last_follow_up > 1460,
  "Alive", "Not yet"
)

AliveAndDead4Years <- AliveAndDead[
  AliveAndDead$AliveStatus == "Alive" | AliveAndDead$deadStatus == "Dead",
]

clinical.w.CIBERSORTx.AliveAndDead <- clinical.w.CIBERSORTx[rownames(AliveAndDead4Years), ]

# Build table of % Alive >4y vs Dead <4y per GES
Percentage.Table <- clinical.w.CIBERSORTx.AliveAndDead %>%
  mutate(
    Status4y = ifelse(
      vital_status == "Dead" & days_to_death < 1460, "Dead < 4 years",
      ifelse(vital_status == "Alive" & days_to_last_follow_up > 1460, "Alive > 4 years", NA)
    )
  ) %>%
  filter(!is.na(Status4y)) %>%
  count(scRNAcluster, Status4y) %>%
  group_by(scRNAcluster) %>%
  mutate(percent = (n / sum(n)) * 100) %>%
  ungroup()

Percentage.Table$scRNAcluster <- factor(
  Percentage.Table$scRNAcluster,
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

# 4. FIGURE 6B – BAR CHART OF 4-YEAR STATUS -----------------------------------
cat("Plotting 4-year alive/dead percentages (Fig 6B)...\n")

p6B <- Percentage.Table %>%
  ggplot(aes(x = scRNAcluster, y = percent, fill = Status4y)) +
  ylab("TCGA-SKCM Samples %") +
  labs(fill = "Vital Status") +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  theme_classic() +
  scale_fill_manual(
    values = c("#F9DC5CFF", "#2460A7FF"),
    breaks = c("Alive > 4 years", "Dead < 4 years"),
    labels = c("Alive > 4 years", "Dead < 4 years")
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18, colour = "black", face = "bold"),
    axis.text.x  = element_text(size = 17, colour = "black", face = "bold", angle = 90),
    axis.text.y  = element_text(size = 13, colour = "black", face = "bold"),
    legend.text  = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 17, face = "bold")
  ) +
  scale_x_discrete(
    labels = c(
      "Anti-apoptosis"                    = "Anti-apoptosis",
      "Immune cell interactions"          = "Immune cells\n interactions",
      "Melanogenesis"                     = "Melanogenesis",
      "Ribosomal biogenesis"              = "Ribosomal\n biogenesis",
      "Extracellular structure organization" = "Extracellular\n structure\n organization",
      "EMT"                               = "EMT"
    )
  )

ggsave("Fig6B_VitalStatus_4Years_Percentage_GESs.png", p6B, width = 8, height = 5, dpi = 1200)

# 5. HELPER: FIT AND PLOT SURVIVAL FOR ONE GES VS REST ------------------------
plot_cluster_survival <- function(cluster_name, cluster_colname, label_cluster, outfile) {
  # cluster_colname: e.g. "Cluster0", "Cluster1", ... if you have such columns
  # here we derive on the fly from scRNAcluster
  
  df <- clinical.w.CIBERSORTx
  
  df[[cluster_colname]] <- ifelse(df$scRNAcluster == cluster_name, cluster_name, "Rest SKCM Patients")
  df[[cluster_colname]] <- factor(df[[cluster_colname]],
                                  levels = c(cluster_name, "Rest SKCM Patients"))
  
  fit_cl <- survfit(
    Surv(overall_survival, deceased) ~ .data[[cluster_colname]],
    data = df
  )
  
  pval_cl <- surv_pvalue(fit_cl, data = df)$pval
  print(paste(cluster_name, "log-rank p =", signif(pval_cl, 3)))
  
  surv_cl <- ggsurvplot(
    fit_cl,
    data       = df,
    pval       = FALSE,
    risk.table = TRUE,
    size       = 2.5,
    legend.labs  = c(label_cluster, "Rest SKCM Patients"),
    legend.title = "",
    font.legend  = c(23, "bold"),
    font.x       = c(25, "bold"),
    font.y       = c(25, "bold"),
    legend       = c(0.59, 0.9),
    palette      = c("blue", "red")
  )
  
  surv_cl$plot <- surv_cl$plot +
    xlab("Time") +
    ylab("Survival probability") +
    annotate(
      "text",
      x      = 2600,
      y      = 0.1,
      label  = paste0("log-rank p = ", signif(pval_cl, 3)),
      cex    = 5,
      vjust  = 0,
      hjust  = 1.1,
      fontface = 4,
      color  = "blue"
    )
  
  ggsave(
    filename = outfile,
    plot     = print(surv_cl$plot, newpage = FALSE),
    device   = "png",
    width    = 6,
    height   = 5,
    dpi      = 1200
  )
}

# 6. FIGURE 6C–H – KM CURVES PER GES -----------------------------------------
cat("Fitting per-GES survival models (Fig 6C–H)...\n")

plot_cluster_survival(
  cluster_name   = "Anti-apoptosis",
  cluster_colname = "Cluster_AA",
  label_cluster  = "Anti-apoptosis",
  outfile        = "Fig6C_AntiApoptosis_vs_Rest.png"
)

plot_cluster_survival(
  cluster_name   = "Immune cell interactions",
  cluster_colname = "Cluster_Immune",
  label_cluster  = "Immune cell interactions",
  outfile        = "Fig6D_ImmuneInteractions_vs_Rest.png"
)

plot_cluster_survival(
  cluster_name   = "Melanogenesis",
  cluster_colname = "Cluster_Melano",
  label_cluster  = "Melanogenesis",
  outfile        = "Fig6E_Melanogenesis_vs_Rest.png"
)

plot_cluster_survival(
  cluster_name   = "Ribosomal biogenesis",
  cluster_colname = "Cluster_Ribo",
  label_cluster  = "Ribosomal biogenesis",
  outfile        = "Fig6F_RibosomalBiogenesis_vs_Rest.png"
)

plot_cluster_survival(
  cluster_name   = "Extracellular structure organization",
  cluster_colname = "Cluster_Extracellular",
  label_cluster  = "Extracellular structure organization",
  outfile        = "Fig6G_Extracellular_vs_Rest.png"
)

plot_cluster_survival(
  cluster_name   = "EMT",
  cluster_colname = "Cluster_EMT",
  label_cluster  = "EMT",
  outfile        = "Fig6H_EMT_vs_Rest.png"
)

cat("FIGURE 6 SCRIPT FINISHED.\n")
cat("Generated files:\n")
cat("- Fig6A_GESs_Global_Survival.png\n")
cat("- Fig6B_VitalStatus_4Years_Percentage_GESs.png\n")
cat("- Fig6C–H_*_vs_Rest.png (one file per GES)\n")
