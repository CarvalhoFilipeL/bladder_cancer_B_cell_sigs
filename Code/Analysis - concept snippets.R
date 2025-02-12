#############################################
## Set working directory and load libraries
#############################################
setwd("H:/My Drive/PostCDDP/ANALYSIS/v4_032024")

# Core libraries for data manipulation, plotting, survival analysis, etc.
library(dplyr)
library(data.table)
library(readr)
library(readxl)
library(survival)
library(ggplot2)
library(ggpubr)
library(flextable)
library(tableone)
library(OmicSelector)    # For correlation plot function (if needed)
library(ComplexHeatmap)  # For heatmap plots
library(circlize)
library(RColorBrewer)
library(meta)
library(metafor)
library(maxstat)
library(gtsummary)

# (Other libraries such as SmartEDA, DataExplorer, etc. are loaded in their sections as needed.)

#############################################
## Define custom functions and settings
#############################################
# Custom grid.draw function for ggsurvplot (if used later)
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

#############################################
## Data Import and Preprocessing
#############################################
# Import merged data from CSV file
merged_d <- fread("Files/merged_d.csv") %>% as.data.frame()

# (Further preprocessing steps can be added here, e.g., filtering, type conversion, etc.)
# For example, if you need to convert IDs:
merged_d$ID <- make.names(merged_d$ID)

#############################################
## Save Data (if needed)
#############################################
fwrite(merged_d, "Files/merged_d.csv")

#############################################
## Heatmap Construction Section
#############################################
# Identify columns for immune, ecotype, and consensus clusters 
# (update column selection as needed)
immune_cols    <- which(names(merged_d) == 'B.cells.naive'):which(names(merged_d) == 'NK_cells')
ecotype_cols   <- which(names(merged_d) == 'CE1'):which(names(merged_d) == 'CE10')
consensus_cols <- which(names(merged_d) == 'LumP'):which(names(merged_d) == 'NE.like')

# Extract the relevant data and compute z-scores (scale columns to mean=0, SD=1)
immune_data     <- as.matrix(merged_d[, immune_cols])
ecotype_data    <- as.matrix(merged_d[, ecotype_cols])
consensus_data  <- as.matrix(merged_d[, consensus_cols])
immune_data_z   <- scale(immune_data)
ecotype_data_z  <- scale(ecotype_data)
consensus_data_z<- scale(consensus_data)

# Replace underscores, dashes, etc. in column names with spaces (for aesthetics)
colnames(immune_data_z)    <- gsub("[._-]", " ", colnames(immune_data_z))
colnames(ecotype_data_z)   <- gsub("[._-]", " ", colnames(ecotype_data_z))
colnames(consensus_data_z) <- gsub("[._-]", " ", colnames(consensus_data_z))

# Define annotation color palettes (customize as needed)
col_Cohort  <- c("GREEK" = "blue", "HM" = "green", "TABER" = "red")
col_T       <- c("T2 or T3" = "yellow", "T4" = "purple")
col_N       <- c("Neg" = "orange", "Pos" = "brown")
col_M       <- c("Neg" = "cyan", "Pos" = "magenta")
col_Platin  <- c("Cisplatin" = "darkblue", "Carboplatin" = "darkred")

# Generate color palettes for other categorical annotations using a helper function
generate_colors <- function(values, palette_name) {
  unique_values <- unique(values[!is.na(values)])
  setNames(colorRampPalette(brewer.pal(8, palette_name))(length(unique_values)), unique_values)
}
col_Carcinoma_Ecotype <- generate_colors(merged_d$`Carcinoma Ecotype`, "Set2")
col_Baylor_subtype    <- generate_colors(merged_d$Baylor.subtype, "Set3")
col_UNC_subtype       <- generate_colors(merged_d$UNC.subtype, "Dark2")
col_CIT_subtype       <- generate_colors(merged_d$CIT.subtype, "Accent")
col_Lund_subtype      <- generate_colors(merged_d$Lund.subtype, "Pastel1")
col_MDA_subtype       <- generate_colors(merged_d$MDA.subtype, "Set1")
col_TCGA_subtype      <- generate_colors(merged_d$TCGA.subtype, "Paired")
col_consensusClass    <- generate_colors(merged_d$consensusClass, "Spectral")

# Create row annotations for the heatmap
row_ha <- rowAnnotation(
  Cohort = merged_d$Cohort,
  Platin = merged_d$Platin,
  Age    = merged_d$Age_at_diagnosis,
  T      = merged_d$T,
  N      = merged_d$N,
  M      = merged_d$M,
  `Carcinoma Ecotype` = merged_d$`Carcinoma Ecotype`,
  `Baylor subtype` = merged_d$Baylor.subtype,
  `UNC subtype` = merged_d$UNC.subtype,
  `CIT subtype` = merged_d$CIT.subtype,
  `Lund subtype` = merged_d$Lund.subtype,
  `MDA subtype` = merged_d$MDA.subtype,
  `TCGA subtype` = merged_d$TCGA.subtype,
  `Consensus Class` = merged_d$consensusClass,
  `OS` = anno_barplot(merged_d$OS_TIME/12, gp = gpar(fill = col_OS <- c("0" = "green", "1" = "red")))
  ,
  col = list(
    Cohort = col_Cohort,
    T = col_T,
    N = col_N,
    M = col_M,
    `Carcinoma Ecotype` = col_Carcinoma_Ecotype,
    `Baylor subtype` = col_Baylor_subtype,
    `UNC subtype` = col_UNC_subtype,
    `CIT subtype` = col_CIT_subtype,
    `Lund subtype` = col_Lund_subtype,
    `MDA subtype` = col_MDA_subtype,
    `TCGA subtype` = col_TCGA_subtype,
    `Consensus Class` = col_consensusClass,
    Platin = col_Platin
  )
)

# Define color function for z-scored heatmap data
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Construct heatmaps for immune infiltration, ecotypes, and consensus clusters
immune_ht <- Heatmap(
  immune_data_z,
  name = "Immune Infiltration",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  left_annotation = row_ha,
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 10)
)
# Draw immune heatmap and record row order for consistency between heatmaps
draw(immune_ht)
row_order <- row_order(immune_ht)

ecotype_ht <- Heatmap(
  ecotype_data_z,
  name = "Ecotypes",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  row_order = row_order,
  column_names_gp = gpar(fontsize = 10)
)

consensus_ht <- Heatmap(
  consensus_data_z,
  name = "Consensus Clusters",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  row_order = row_order,
  column_names_gp = gpar(fontsize = 10)
)

# Combine and draw the heatmaps, saving output as PNG
ht_list <- immune_ht + ecotype_ht + consensus_ht
png("Plots/Heatmap.png", width = 20, height = 12, units = "in", res = 300)
draw(ht_list, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

#############################################
## Exploratory Data Analysis (EDA)
#############################################
# Using dataMaid, SmartEDA, and DataExplorer to generate reports

# Uncomment and install if packages are not yet installed:
# install.packages("dataMaid")
library(dataMaid)
makeDataReport(merged_d, output = "html", replace = TRUE, file = "EDA/dataMaid.html")

# install.packages("SmartEDA")
library(SmartEDA)
ExpReport(merged_d, op_file = 'EDA/smartEDA.html')

# install.packages("DataExplorer")
library(DataExplorer)
create_report(merged_d, output_file = "EDA/DataExplorer.html")

#############################################
## Table 1 Generation using tableone package
#############################################
temp_table <- merged_d %>% select(-Cisplatin, -Carboplatin)
temp_table$Age_at_diagnosis <- as.numeric(temp_table$Age_at_diagnosis)
temptable <- CreateTableOne(
  data = temp_table,
  vars = c("Age_at_diagnosis", "StageIV", "T", "N", "M", "MVAC", "Platin"),
  includeNA = TRUE,
  strata = "Cohort"
)
temptable_mat <- print(temptable, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(temptable_mat, file = "Files/table1.csv")
flextable(as.data.frame(temptable_mat)) %>% save_as_docx(path = "Files/Table1.docx")

#############################################
## Survival Analyses: Descriptive and OS Comparisons
#############################################
# Example: Kaplan-Meier survival curves comparing cohorts
s1 <- survfit(Surv(OS_TIME, OS) ~ Cohort, data = merged_d)
temphtml <- tbl_survfit(s1, times = c(12, 24, 60), label_header = "**{time}-month survival**")
temphtml %>% as_flex_table() %>% save_as_docx(path = "Files/OStimes_Survival_vs_Cohort.docx")

# Save KM plot with risk table
km_plot <- ggsurvplot(s1, pval = TRUE, risk.table = TRUE, break.time.by = 12, xlim = c(0,84))
ggsave(filename = "Plots/OS_vs_cohorts.png", plot = km_plot, width = 8, height = 6, dpi = 300)

# Univariable Cox model for Cohort effect
cox_model <- coxph(Surv(OS_TIME, OS) ~ Cohort, data = merged_d)
tbl_regression(cox_model, exponentiate = TRUE) %>% add_global_p() %>%
  as_flex_table() %>% save_as_docx(path = "Files/UnivariableOS_Cohort.docx")

# Similar analyses for other variables (e.g., T, N, M, StageIV, Platin, etc.) follow here...
# (Repeat similar steps: create survfit, generate KM plot, run coxph, save tables and plots.)

#############################################
## Analysis with Continuous Variables
#############################################
# Example for Age_at_diagnosis using maxstat for optimal cutpoint
tempsplit <- maxstat.test(Surv(OS_TIME, OS) ~ Age_at_diagnosis, data = merged_d,
                          smethod = "LogRank", pmethod = "Lau94")
merged_d$Group <- as.factor(
  ifelse(merged_d$Age_at_diagnosis > tempsplit$estimate,
         paste0("Age at diagnosis >", round(tempsplit$estimate, 0)),
         paste0("Age at diagnosis <=", round(tempsplit$estimate, 0)))
)
# Survival analysis with new group variable
s1 <- survfit(Surv(OS_TIME, OS) ~ Group, data = merged_d)
s1 %>% tbl_survfit(times = c(12, 24, 60), label_header = "**{time}-month survival**") %>%
  as_flex_table() %>% save_as_docx(path = "Files/OStime_Survival_vs_AgeGroups.docx")
km_age <- ggsurvplot(s1, pval = TRUE, risk.table = TRUE, break.time.by = 12, xlim = c(0,84))
ggsave("Plots/OS_vs_AgeGroups.png", plot = km_age, width = 8, height = 6, dpi = 300)
m_age <- coxph(Surv(OS_TIME, OS) ~ Age_at_diagnosis, data = merged_d)
tbl_regression(m_age, exponentiate = TRUE) %>% as_flex_table() %>%
  save_as_docx(path = "Files/UnivariableOS_Age_at_diagnosis.docx")

#############################################
## Metaanalysis Section
#############################################
# The script below iterates over immune parameters and cohorts to run Cox models,
# extract HRs and confidence intervals, then perform meta-analyses.
# Adjust column selection (immuoparameters) as required.
immuoparameters <- colnames(merged_d)[c(immune_cols, ecotype_cols, consensus_cols)]
if (!dir.exists("reports")) dir.create("reports")
if (!dir.exists("forests")) dir.create("forests")

# Ensure Platin is a factor and check data structure
merged_d$Platin <- as.factor(merged_d$Platin)

# Loop over cohorts and immune parameters to build univariable Cox models
kohorty <- unique(as.character(merged_d$Cohort))
adjustedHRs <- data.frame()

for (kohorta in kohorty) {
  cohort_data <- merged_d[merged_d$Cohort == kohorta, ]
  wynikidf <- data.frame()
  for (selected in immuoparameters) {
    try({
      # Select variables and remove missing cases
      tempxx <- cohort_data %>% select(OS_TIME, OS, Age_at_diagnosis, T, M, Platin, all_of(selected))
      tempxx <- tempxx[complete.cases(tempxx), ]
      
      # Build formula and run Cox model
      fmla <- as.formula(paste("Surv(OS_TIME, OS) ~ Age_at_diagnosis + T + M +", selected))
      m <- coxph(fmla, data = tempxx)
      
      # Extract model summary statistics
      wynikidf[selected, "parameter"] <- selected
      wynikidf[selected, "loglik"] <- logLik(m)
      wynikidf[selected, "HR"] <- summary(m)$conf.int[4, 1]
      wynikidf[selected, "LowerCI"] <- summary(m)$conf.int[4, 3]
      wynikidf[selected, "UpperCI"] <- summary(m)$conf.int[4, 4]
      wynikidf[selected, "cindex"] <- summary(m)$concordance["C"]
      wynikidf[selected, "pval"] <- summary(m)$coefficients[4, 5]
      wynikidf[selected, "cohort"] <- kohorta
      wynikidf[selected, "n"] <- nrow(cohort_data)
    })
  }
  adjustedHRs <- rbind(adjustedHRs, wynikidf)
}
fwrite(adjustedHRs, "adjustedHRs.csv")

# Meta-analysis for each immune parameter
wyniki_metaanaliz <- data.frame()
for (parametr in immuoparameters) {
  temp <- adjustedHRs %>% filter(parameter == parametr)
  temp$logHR <- log(temp$HR)
  temp$logUB <- log(temp$UpperCI)
  temp$logLB <- log(temp$LowerCI)
  temp$selogHR <- (temp$logUB - temp$logLB) / (2 * 1.96)
  
  metaanalysis <- metagen(temp$logHR, temp$selogHR, studlab = temp$cohort, data = temp, sm = "HR")
  if (!is.null(metaanalysis)) {
    sink(paste0("reports/", parametr, ".txt"))
    print(summary(metaanalysis))
    sink()
    
    wyniki_metaanaliz[parametr, "parameter"] <- parametr
    wyniki_metaanaliz[parametr, "HR_fixed"] <- exp(metaanalysis$TE.fixed)
    wyniki_metaanaliz[parametr, "LowerCI_fixed"] <- exp(metaanalysis$lower.fixed)
    wyniki_metaanaliz[parametr, "UpperCI_fixed"] <- exp(metaanalysis$upper.fixed)
    wyniki_metaanaliz[parametr, "p_fixed"] <- metaanalysis$pval.fixed
    wyniki_metaanaliz[parametr, "I2"] <- metaanalysis$I2
    wyniki_metaanaliz[parametr, "HR_random"] <- exp(metaanalysis$TE.random)
    wyniki_metaanaliz[parametr, "LowerCI_random"] <- exp(metaanalysis$lower.random)
    wyniki_metaanaliz[parametr, "UpperCI_random"] <- exp(metaanalysis$upper.random)
    wyniki_metaanaliz[parametr, "p_random"] <- metaanalysis$pval.random
    
    pdf(paste0("forests/", parametr, ".pdf"), width = 20)
    forest(metaanalysis, title = parametr)
    dev.off()
    
    try({
      a <- metabias(metaanalysis, method = "rank", k.min = 2)
      wyniki_metaanaliz[parametr, "BeggMazumdar_p"] <- a$pval
    })
    try({
      a <- metabias(metaanalysis, method = "linreg", k.min = 2)
      wyniki_metaanaliz[parametr, "Egger_p"] <- a$pval
    })
  }
}
fwrite(wyniki_metaanaliz, "metaanalyses.csv")

#############################################
## Custom Plots (Boxplots and Correlations)
#############################################
# Define a custom color palette for cohorts
custom_colors <- c("GREEK" = "#F8766D", "HM" = "#00BA38", "TABER" = "#619CFF")

# Boxplot: TILs by Cohort
png("Plots/Boxplot_Cohort_TILs.png", width = 6, height = 4, units = "in", res = 300)
ggboxplot(merged_d, x = "Cohort", y = "TILs", add = "jitter", fill = "Cohort") +
  scale_fill_manual(values = custom_colors) +
  geom_pwc(tip.length = 0, label = "p = {p.adj.format}", method = "t_test") +
  theme_minimal()
dev.off()

# Boxplot: B cells by Cohort
png("Plots/Boxplot_Cohort_Bcells.png", width = 6, height = 4, units = "in", res = 300)
ggboxplot(merged_d, x = "Cohort", y = "B_cells", add = "jitter", fill = "Cohort") +
  scale_fill_manual(values = custom_colors) +
  geom_pwc(tip.length = 0, label = "p = {p.adj.format}", method = "t_test") +
  ylab("B cells") +
  theme_minimal()
dev.off()

# Boxplot: Lymphoid lineage by Cohort
png("Plots/Boxplot_Cohort_Lymphoid.png", width = 6, height = 4, units = "in", res = 300)
ggboxplot(merged_d, x = "Cohort", y = "Lymphoid", add = "jitter", fill = "Cohort") +
  scale_fill_manual(values = custom_colors) +
  geom_pwc(tip.length = 0, label = "p = {p.adj.format}", method = "t_test") +
  ylab("Lymphoid lineage") +
  theme_minimal()
dev.off()

# Boxplot: B cells memory by Cohort
png("Plots/Boxplot_Cohort_BcellsMem.png", width = 6, height = 4, units = "in", res = 300)
ggboxplot(merged_d, x = "Cohort", y = "B.cells.memory", add = "jitter", fill = "Cohort") +
  scale_fill_manual(values = custom_colors) +
  geom_pwc(tip.length = 0, label = "p = {p.adj.format}", method = "t_test") +
  ylab("B cells memory") +
  theme_minimal()
dev.off()

# Correlation plot using OmicSelector (if available)
png("Plots/Correlation_CE10_BcellsMem.png", width = 6, height = 5, units = "in", res = 300)
OmicSelector_correlation_plot(merged_d$CE10, merged_d$B.cells.memory,
                              "CE10 score", "B cells memory", "", yx = FALSE)
dev.off()

# Additional custom correlation plots (CE6, CE7, LumNS, Neutrophils, etc.)
png("Plots/Correlation_CE6_BcellsMem.png", width = 6, height = 5, units = "in", res = 300)
OmicSelector_correlation_plot(merged_d$CE6, merged_d$B.cells.memory,
                              "CE6 score", "B cells memory", "", yx = FALSE)
dev.off()

png("Plots/Correlation_CE7_BcellsMem.png", width = 6, height = 5, units = "in", res = 300)
OmicSelector_correlation_plot(merged_d$CE7, merged_d$B.cells.memory,
                              "CE7 score", "B cells memory", "", yx = FALSE)
dev.off()

png("Plots/Correlation_LumNS_BcellsMem.png", width = 6, height = 5, units = "in", res = 300)
OmicSelector_correlation_plot(merged_d$LumNS, merged_d$B.cells.memory,
                              "LumNS", "B cells memory", "", yx = FALSE, metoda = "spearman")
dev.off()

png("Plots/Correlation_Neutrophils_BcellsMem.png", width = 6, height = 5, units = "in", res = 300)
OmicSelector_correlation_plot(merged_d$Neutrophils, merged_d$B.cells.memory,
                              "Neutrophils", "B cells memory", "", yx = FALSE)
dev.off()

#############################################
## Final Plot: Scatter Plot Example
#############################################
# Example of a scatter plot with linear regression line and Spearman correlation
# Adjust axis limits and labels as needed.
plot1 <- ggplot(merged_d, aes(x = bcell_signature1, y = B.cells.memory)) +
  geom_point() +
  geom_smooth(method = "lm", color = "darkred", fill = "gray") +
  stat_cor(method = "spearman", label.x = 0, label.y = 1, p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(title = "", x = "IFNG Expression / Total CD8 Fraction", y = "B-cell memory (CIBERSORTx)") +
  ylim(0, 1) +
  theme_bw()
plot1

# If additional plots are needed, combine them using ggpubr::ggarrange
