setwd("H:/My Drive/PostCDDP/ANALYSIS/v4_032024")

library(dplyr)
library(flextable)
# rspm::enable()

# add method to grid.draw
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

merged_d = data.table::fread("../v3/data.csv") %>% as.data.frame()
merged_d_full = data.table::fread("../v3/data.csv")%>% as.data.frame()
merged_d = merged_d %>% dplyr::filter(Cohort == "GREEK" | Cohort == "HM" | Cohort == "TABER") %>% dplyr::filter(Platin != "None")
merged_d$Platin = as.character(merged_d$Platin)
merged_d$Cohort = as.character(merged_d$Cohort)

merged_d = merged_d %>% mutate_if(is.character, as.factor)

table(merged_d$T == "T2 or T3" & merged_d$N == "Neg" & merged_d$M == "Neg")

merged_d$Operable = merged_d$T == "T2 or T3" & merged_d$N == "Neg" & merged_d$M == "Neg"
table(merged_d$Operable, merged_d$Cohort)

library(readxl)
X41467_2020_18640_MOESM4_ESM_1_ <- read_excel("TABERcheck/41467_2020_18640_MOESM4_ESM (1).xlsx", 
                                              sheet = "data")
got_neoadjuvant_ids = X41467_2020_18640_MOESM4_ESM_1_$ID[X41467_2020_18640_MOESM4_ESM_1_$Neoadjucant == "Yes"]

merged_d$Neoadjuvant = ifelse(merged_d$ID %in% got_neoadjuvant_ids, "Yes", "No")
table(merged_d$Neoadjuvant)

merged_d = merged_d %>% dplyr::filter(Neoadjuvant != "Yes") %>% dplyr::select(-Neoadjuvant)
data.table::fwrite(merged_d, "Files/merged_d.csv")

# Derived
merged_d$TILs = merged_d$B.cells.naive + merged_d$B.cells.memory + merged_d$Plasma.cells + merged_d$T.cells.CD8 + merged_d$T.cells.CD4.naive + merged_d$T.cells.CD4.memory.resting +
  merged_d$T.cells.CD4.memory.activated + merged_d$T.cells.follicular.helper + merged_d$T.cells.regulatory..Tregs. + merged_d$T.cells.gamma.delta +
  merged_d$NK.cells.resting + merged_d$NK.cells.activated + merged_d$Monocytes + merged_d$Macrophages.M0 +
  merged_d$Macrophages.M1 + merged_d$Macrophages.M2 + merged_d$Dendritic.cells.resting + merged_d$Dendritic.cells.activated + merged_d$Mast.cells.resting + merged_d$Mast.cells.activated
merged_d$Myeloid = merged_d$Eosinophils + merged_d$Monocytes + merged_d$Macrophages.M0 + merged_d$Macrophages.M1 + merged_d$Macrophages.M2 + merged_d$Neutrophils + merged_d$Mast.cells.resting + merged_d$Mast.cells.activated
merged_d$Monocytes_derivatives = merged_d$Macrophages.M0 + merged_d$Macrophages.M1 + merged_d$Macrophages.M2 + merged_d$Dendritic.cells.activated + merged_d$Dendritic.cells.resting + merged_d$Monocytes
merged_d$Macrophages = merged_d$Macrophages.M0 + merged_d$Macrophages.M1 + merged_d$Macrophages.M2


merged_d$Lymphoid = merged_d$B.cells.memory + merged_d$B.cells.naive + merged_d$Plasma.cells + merged_d$T.cells.CD4.memory.activated + merged_d$T.cells.CD4.memory.resting + merged_d$T.cells.CD4.naive + merged_d$T.cells.CD8 + merged_d$T.cells.follicular.helper + merged_d$T.cells.gamma.delta + merged_d$T.cells.regulatory..Tregs. + merged_d$NK.cells.activated + merged_d$NK.cells.resting
merged_d$B_cells = merged_d$B.cells.memory + merged_d$B.cells.naive
merged_d$B_and_plasma_cells = merged_d$B.cells.memory + merged_d$B.cells.naive + merged_d$Plasma.cells
merged_d$Tcells_CD4 = merged_d$T.cells.CD4.memory.activated + merged_d$T.cells.CD4.memory.resting + merged_d$T.cells.CD4.naive
merged_d$Tcells = merged_d$T.cells.CD4.memory.activated + merged_d$T.cells.CD4.memory.resting + merged_d$T.cells.CD4.naive + merged_d$T.cells.CD8 + + merged_d$T.cells.follicular.helper + + merged_d$T.cells.gamma.delta + merged_d$T.cells.regulatory..Tregs.
merged_d$NK_cells = merged_d$NK.cells.activated + merged_d$NK.cells.resting

data.table::fwrite(merged_d, "Files/merged_d.csv")

# EcoTyper ----
temp1 = data.table::fread("Files/ecotyper_HM+GREEK/Carcinoma_Ecotypes/Ecotype_Abundance.txt")
temp1_1 = data.table::fread("Files/ecotyper_HM+GREEK/Carcinoma_Ecotypes/Ecotype_Assignment.txt")

temp_df = dplyr::full_join(temp1, temp1_1, by="ID")

temp2 = data.table::fread("Files/ecotyper_TABER/Carcinoma_Ecotypes/Ecotype_Abundance.txt")
temp2_1 = data.table::fread("Files/ecotyper_TABER/Carcinoma_Ecotypes/Ecotype_Assignment.txt")

temp_df2 = dplyr::full_join(temp2, temp2_1, by="ID")

temp = rbind(temp_df, temp_df2)
merged_d$ID_orginal = merged_d$ID

merged_d$ID = make.names(merged_d$ID)

merged_d = dplyr::left_join(merged_d, temp, by = "ID")
merged_d


table(merged_d$`Carcinoma Ecotype`)
colnames(merged_d)[colnames(merged_d) == "Carcionoma Ecotype"] = "CE"
data.table::fwrite(merged_d, "Files/merged_d.csv")

# Classes ----
temp2 = read.csv("../v1_Dec2022/HM+GREEK/mRNA_Clusters/clusters.csv")
temp3 = read.csv("../v1_Dec2022/HM+GREEK/mRNA_Clusters/subtypes.csv")

temp1 = cbind(temp3, temp2)

temp6 = read.csv("../v1_Dec2022/TABER/mRNA_Clusters/clusters.csv")
temp7 = read.csv("../v1_Dec2022/TABER/mRNA_Clusters/subtypes.csv")

temp4 = cbind(temp7,temp6)

temp8 = rbind(temp1, temp4)
table(temp8$ID)
temp8$ID = make.names(temp8$ID)

merged_d = dplyr::left_join(merged_d, temp8, by = "ID")
merged_d = dplyr::select(merged_d, -cor_pval, -separationLevel, -Operable)
data.table::fwrite(merged_d, "Files/merged_d.csv")

# Heatmap ----
# Install and load necessary packages
# install.packages("ComplexHeatmap")
# install.packages("circlize")
# install.packages("RColorBrewer")

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Identify columns for each heatmap
immune_cols <- which(names(merged_d) == 'B.cells.naive'):which(names(merged_d) == 'NK_cells')
ecotype_cols <- which(names(merged_d) == 'CE1'):which(names(merged_d) == 'CE10')
consensus_cols <- which(names(merged_d) == 'LumP'):which(names(merged_d) == 'NE.like')

# Extract data for each heatmap
immune_data <- as.matrix(merged_d[, immune_cols])
ecotype_data <- as.matrix(merged_d[, ecotype_cols])
consensus_data <- as.matrix(merged_d[, consensus_cols])

# Z-score the data (scale each column to have mean 0 and SD 1)
immune_data_z <- scale(immune_data)
ecotype_data_z <- scale(ecotype_data)
consensus_data_z <- scale(consensus_data)

# Replace underscores and dashes with spaces in column names
colnames(immune_data_z) <- gsub("[._-]", " ", colnames(immune_data_z))
colnames(ecotype_data_z) <- gsub("[._-]", " ", colnames(ecotype_data_z))
colnames(consensus_data_z) <- gsub("[._-]", " ", colnames(consensus_data_z))

# Also replace in the row annotations' labels if necessary
# Create a named vector to map old names to new names with spaces
annotation_labels <- c(
  "Age_at_diagnosis" = "Age at diagnosis",
  "OS_TIME" = "OS TIME",
  "OS" = "OS",
  "StageIV" = "Stage IV",
  "Carcinoma Ecotype" = "Carcinoma Ecotype",
  "Baylor_subtype" = "Baylor subtype",
  "UNC_subtype" = "UNC subtype",
  "CIT_subtype" = "CIT subtype",
  "Lund_subtype" = "Lund subtype",
  "MDA_subtype" = "MDA subtype",
  "TCGA_subtype" = "TCGA subtype",
  "consensusClass" = "Consensus Class"
)

# Define colors for categorical annotations
col_Cohort <- c("GREEK" = "blue", "HM" = "green", "TABER" = "red")
col_T <- c("T2 or T3" = "yellow", "T4" = "purple")
col_N <- c("Neg" = "orange", "Pos" = "brown")
col_M <- c("Neg" = "cyan", "Pos" = "magenta")
col_StageIV <- c("No" = "grey", "Yes" = "black")
col_Platin = c("Cisplatin" = "darkblue", "Carboplatin" = "darkred")

# Define colors for OS
col_OS <- c("0" = "green", "1" = "red")  # Assuming 0 and 1 represent survival status
OS_char <- as.character(merged_d$OS)
colors_OS <- col_OS[OS_char]

# Function to generate color palettes for categorical variables
generate_colors <- function(values, palette_name) {
  unique_values <- unique(values)
  unique_values <- unique_values[!is.na(unique_values)]  # Remove NA values
  colors <- setNames(colorRampPalette(brewer.pal(8, palette_name))(length(unique_values)), unique_values)
  return(colors)
}

# Generate color mappings for new categorical variables
col_Carcinoma_Ecotype <- generate_colors(merged_d$`Carcinoma Ecotype`, "Set2")
col_Baylor_subtype <- generate_colors(merged_d$Baylor.subtype, "Set3")
col_UNC_subtype <- generate_colors(merged_d$UNC.subtype, "Dark2")
col_CIT_subtype <- generate_colors(merged_d$CIT.subtype, "Accent")
col_Lund_subtype <- generate_colors(merged_d$Lund.subtype, "Pastel1")
col_MDA_subtype <- generate_colors(merged_d$MDA.subtype, "Set1")
col_TCGA_subtype <- generate_colors(merged_d$TCGA.subtype, "Paired")
col_consensusClass <- generate_colors(merged_d$consensusClass, "Spectral")

# Create row annotations
row_ha <- rowAnnotation(
  Cohort = merged_d$Cohort,
  Platin = merged_d$Platin,
  
  
  `Age` = merged_d$Age_at_diagnosis,
  T = merged_d$T,
  N = merged_d$N,
  M = merged_d$M,
  `Carcinoma Ecotype` = merged_d$`Carcinoma Ecotype`,
  `Baylor subtype` = merged_d$Baylor.subtype,
  `UNC subtype` = merged_d$UNC.subtype,
  `CIT subtype` = merged_d$CIT.subtype,
  `Lund subtype` = merged_d$Lund.subtype,
  `MDA subtype` = merged_d$MDA.subtype,
  `TCGA subtype` = merged_d$TCGA.subtype,
  `Consensus Class` = merged_d$consensusClass,
  `OS` = anno_barplot(merged_d$OS_TIME/12, 
                      gp = gpar(fill = colors_OS)),
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

# Define the color function for z-scored data
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Construct the immune infiltration heatmap with clustering using z-scored data
immune_ht <- Heatmap(immune_data_z,
                     name = "Immune Infiltration",
                     col = col_fun,
                     cluster_rows = TRUE,
                     cluster_columns = TRUE,
                     left_annotation = row_ha,
                     show_row_names = FALSE,
                     column_names_gp = gpar(fontsize = 10))

# Draw the immune heatmap to get the row order
draw(immune_ht)
row_order <- row_order(immune_ht)

# Create ecotype heatmap using z-scored data and same row order
ecotype_ht <- Heatmap(ecotype_data_z,
                      name = "Ecotypes",
                      col = col_fun,
                      cluster_rows = FALSE,
                      cluster_columns = TRUE,
                      show_row_names = FALSE,
                      row_order = row_order,
                      column_names_gp = gpar(fontsize = 10))

# Create consensus cluster heatmap using z-scored data and same row order
consensus_ht <- Heatmap(consensus_data_z,
                        name = "Consensus Clusters",
                        col = col_fun,
                        cluster_rows = FALSE,
                        cluster_columns = TRUE,
                        show_row_names = FALSE,
                        row_order = row_order,
                        column_names_gp = gpar(fontsize = 10))

# Combine the heatmaps
ht_list <- immune_ht + ecotype_ht + consensus_ht

# Draw the combined heatmap
png("Plots/Heatmap.png", width = 20, height = 12, units = "in", res = 300)
draw(ht_list, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()




# EDA ----
temp = merged_d

install.packages("dataMaid")
library(dataMaid)
makeDataReport(temp, output = "html", replace = TRUE, file = "EDA/dataMaid.html")

install.packages("SmartEDA")
library(SmartEDA)
ExpReport(temp,op_file='EDA/smartEDA.html')

install.packages("DataExplorer")
library(DataExplorer)
DataExplorer::create_report(temp, output_file = "EDA/DataExplorer.html")

# Table 1 ----
temp = merged_d
temp = dplyr::select(temp, -Cisplatin, -Carboplatin)
colnames(temp)
str(temp)

library(tableone)
temp$Age_at_diagnosis = as.numeric(temp$Age_at_diagnosis)
temptable = CreateTableOne(data = temp, vars = c("Age_at_diagnosis","StageIV", "T", "N", "M", "MVAC", "Platin"), includeNA = T, strata = "Cohort")
temptable_mat <- print(temptable, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(temptable_mat, file = "Files/table1.csv")
temptable_mat
flextable(as.data.frame(temptable_mat)) %>% flextable::save_as_docx(path="Files/Table1.docx")

# Descriptions ----
# Load necessary libraries
library(survival)
library(dplyr)

temps = dplyr::filter(temp, Cohort == "TABER")
median(temps$OS_TIME)
max(temps$OS_TIME)
table(temps$OS)
survfit(Surv(OS_TIME, OS) ~ 1, data = temps)
summary(survfit(Surv(OS_TIME, OS) ~ 1, data = temps))



# OS comparisons ----
library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(ggpubr)
library(survminer)
library(dplyr)
compare_survival_at_time = function(fit, time){
  threeYr = summary(fit,times=time)
  #difference in survival at 3 years between edema=0 and edemo=1 (for example) is
  threeYr$surv[1] - threeYr$surv[3]
  #the standard error of this is
  diffSE <- sqrt(threeYr$std.err[3]^2 + threeYr$std.err[1]^2)
  #a 95% CI for the diff is
  threeYr$surv[1] - threeYr$surv[3] - 1.96 *diffSE
  threeYr$surv[1] - threeYr$surv[3] + 1.96 *diffSE
  #a z-test test statistic is
  zStat <- (threeYr$surv[1] - threeYr$surv[3])/diffSE
  #and a two-sided p-value testing that the diff. is 0 is
  return(2*pnorm(abs(zStat), lower.tail=FALSE))
}

#
s1 <- survfit(Surv(OS_TIME, OS) ~ Cohort, data = temp)
s1

temphtml = s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" )
temphtml
temphtml %>% as_flex_table() %>% flextable::save_as_docx(path="Files/OStimes_Survival_vs_Cohort.docx")

a = ggsurvplot(s1, pval = T, risk.table = T, break.time.by=12, xlim = c(0,84))
k = ggsurvplot(s1, pval = T, risk.table = F, break.time.by=12, xlim = c(0,84))
ggsave(file = "Plots/OS_vs_cohorts.png", plot= a, width = 8, height = 6, dpi = 300)

m = coxph(Surv(OS_TIME, OS) ~ Cohort, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% add_global_p()
tbl_regression(m, exponentiate = T) %>% add_global_p() %>% as_flex_table() %>% flextable::save_as_docx(path="Files/UnivariableOS_Cohort.docx")


# Types
#
temp$CE = temp$`Carcinoma Ecotype`
s1 <- survfit(Surv(OS_TIME, OS) ~ CE, data = temp)
s1

temphtml = s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" )
temphtml
temphtml %>% as_flex_table() %>% flextable::save_as_docx(path="Files/OStimes_Survival_vs_CE.docx")

a = ggsurvplot(s1, pval = T, risk.table = T, break.time.by=12, xlim = c(0,84))
k = ggsurvplot(s1, pval = T, risk.table = F, break.time.by=12, xlim = c(0,84))
ggsave(file = "Plots/OS_vs_CE.png", plot= a, width = 8, height = 8, dpi = 300)

m = coxph(Surv(OS_TIME, OS) ~ CE, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% add_global_p()
tbl_regression(m, exponentiate = T) %>% add_global_p() %>% as_flex_table() %>% flextable::save_as_docx(path="Files/UnivariableOS_CE.docx")

#
s1 <- survfit(Surv(OS_TIME, OS) ~ Platin, data = temp)
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" )
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" ) %>% as_flex_table() %>% flextable::save_as_docx(path="Files/OStime_Survival_vs_Platin.docx")
e = ggsurvplot(s1, pval = T, risk.table = T, break.time.by=12, xlim = c(0,84))
e

ggsave(file = "Plots/OS_vs_Platin.png", plot= e, width = 8, height = 6, dpi = 300)

m = coxph(Surv(OS_TIME, OS) ~ Platin, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% add_global_p()
tbl_regression(m, exponentiate = T) %>% add_global_p() %>% as_flex_table() %>% flextable::save_as_docx(path="Files/UnivariableOS_Platin.docx")

## T ----
s1 <- survfit(Surv(OS_TIME, OS) ~ T, data = temp)
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" )
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" ) %>% as_flex_table() %>% flextable::save_as_docx(path="Files/OStime_Survival_vs_T.docx")
e = ggsurvplot(s1, pval = T, risk.table = T, break.time.by=12, xlim = c(0,84))
e

ggsave(file = "Plots/OS_vs_T.png", plot= e, width = 8, height = 6, dpi = 300)

m = coxph(Surv(OS_TIME, OS) ~ T, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% add_global_p()
tbl_regression(m, exponentiate = T) %>% add_global_p() %>% as_flex_table() %>% flextable::save_as_docx(path="Files/UnivariableOS_T.docx")

## N ----
temp$N
s1 <- survfit(Surv(OS_TIME, OS) ~ N, data = temp)
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" )
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" ) %>% as_flex_table() %>% flextable::save_as_docx(path="Files/OStime_Survival_vs_N.docx")
e = ggsurvplot(s1, pval = T, risk.table = T, break.time.by=12, xlim = c(0,84))
e

ggsave(file = "Plots/OS_vs_N.png", plot= e, width = 8, height = 6, dpi = 300)

m = coxph(Surv(OS_TIME, OS) ~ N, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% add_global_p()
tbl_regression(m, exponentiate = T) %>% add_global_p() %>% as_flex_table() %>% flextable::save_as_docx(path="Files/UnivariableOS_N.docx")

## M ----
temp$M
s1 <- survfit(Surv(OS_TIME, OS) ~ M, data = temp)
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" )
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" ) %>% as_flex_table() %>% flextable::save_as_docx(path="Files/OStime_Survival_vs_M.docx")
e = ggsurvplot(s1, pval = T, risk.table = T, break.time.by=12, xlim = c(0,84))
e

ggsave(file = "Plots/OS_vs_M.png", plot= e, width = 8, height = 6, dpi = 300)

m = coxph(Surv(OS_TIME, OS) ~ M, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% add_global_p()
tbl_regression(m, exponentiate = T) %>% add_global_p() %>% as_flex_table() %>% flextable::save_as_docx(path="Files/UnivariableOS_M.docx")

## StageIV ----
temp$StageIV
s1 <- survfit(Surv(OS_TIME, OS) ~ StageIV, data = temp)
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" )
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" ) %>% as_flex_table() %>% flextable::save_as_docx(path="Files/OStime_Survival_vs_StageIV.docx")
e = ggsurvplot(s1, pval = T, risk.table = T, break.time.by=12, xlim = c(0,84))
e

ggsave(file = "Plots/OS_vs_StageIV.png", plot= e, width = 8, height = 6, dpi = 300)

m = coxph(Surv(OS_TIME, OS) ~ StageIV, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% add_global_p()
tbl_regression(m, exponentiate = T) %>% add_global_p() %>% as_flex_table() %>% flextable::save_as_docx(path="Files/UnivariableOS_StageIV.docx")


## Platin ----
temp$Platin
s1 <- survfit(Surv(OS_TIME, OS) ~ Platin, data = temp)
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" )
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" ) %>% as_flex_table() %>% flextable::save_as_docx(path="Files/OStime_Survival_vs_Platin.docx")
e = ggsurvplot(s1, pval = T, risk.table = T, break.time.by=12, xlim = c(0,84))
e

ggsave(file = "Plots/OS_vs_Platin.png", plot= e, width = 8, height = 6, dpi = 300)

m = coxph(Surv(OS_TIME, OS) ~ Platin, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% add_global_p()
tbl_regression(m, exponentiate = T) %>% add_global_p() %>% as_flex_table() %>% flextable::save_as_docx(path="Files/UnivariableOS_Platin.docx")


## MVAC ----
temp$MVAC
s1 <- survfit(Surv(OS_TIME, OS) ~ MVAC, data = temp)
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" )
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" ) %>% as_flex_table() %>% flextable::save_as_docx(path="Files/OStime_Survival_vs_MVAC.docx")
e = ggsurvplot(s1, pval = T, risk.table = T, break.time.by=12, xlim = c(0,84))
e

ggsave(file = "Plots/OS_vs_MVAC.png", plot= e, width = 8, height = 6, dpi = 300)

m = coxph(Surv(OS_TIME, OS) ~ MVAC, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% add_global_p()
tbl_regression(m, exponentiate = T) %>% add_global_p() %>% as_flex_table() %>% flextable::save_as_docx(path="Files/UnivariableOS_MVAC.docx")

## Continous ----
nazwa = "Age at diagnosis"
var_name = "Age_at_diagnosis"

library(maxstat)
tempsplit = maxstat.test(Surv(OS_TIME, OS) ~ Age_at_diagnosis,  data=temp, smethod="LogRank", pmethod="Lau94")

temp$Group = as.factor(ifelse(temp$Age_at_diagnosis>tempsplit$estimate,paste0(nazwa,">",round(tempsplit$estimate,0)),paste0(nazwa,"<=",round(tempsplit$estimate,0))))
table(temp$Group)

s1 <- survfit(Surv(OS_TIME, OS) ~ Group, data = temp)
s1
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-month survival**" )
s1 %>% tbl_survfit(times = c(12,24,60), label_header = "**{time}-year survival**") %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0("Files/OStime_Survival_vs_",nazwa,".docx"))

f = ggsurvplot(s1, pval = T, risk.table = T, break.time.by=12, xlim = c(0,84))
f

ggsave(file = "Plots/OS_vs_AgeGroups.png", plot= f, width = 8, height = 6, dpi = 300)

m = coxph(Surv(OS_TIME, OS)  ~ Age_at_diagnosis, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0("Files/UnivariableOS_",nazwa,".docx"))
tbl_regression(m, exponentiate = T)


### CE10 ----
nazwa = "CE10"
var_name = "CE10"

temp = merged_d
library(maxstat)
tempsplit = maxstat.test(Surv(OS_TIME, OS) ~ CE10,  data=temp, smethod="LogRank", pmethod="Lau94")

temp$Group = as.factor(ifelse(temp$CE10>tempsplit$estimate,paste0(nazwa,">",round(tempsplit$estimate,2)),paste0(nazwa,"<=",round(tempsplit$estimate,2))))

s1 <- survfit(Surv(OS_TIME, OS) ~ Group, data = temp)
s1 %>% tbl_survfit(times = c(1,2,5), label_header = "**{time}-year survival**") %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0("Files/OStime_Survival_vs_",nazwa,".docx"))
png(paste0("Plots/OS_vs_",nazwa,".png"), width = 2000, height = 1500, res = 300)
ggsurvplot(s1, pval = T, risk.table = T, legend.labs = levels(as.factor(temp$Group)), break.time.by = 12, xlim=c(0,60))
dev.off()
m = coxph(Surv(OS_TIME, OS)  ~ Group, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0("Files/UnivariableOS_",nazwa,".docx"))


nazwa = "CE9"
var_name = "CE9"

temp = merged_d
library(maxstat)
tempsplit = maxstat.test(Surv(OS_TIME, OS) ~ CE9,  data=temp, smethod="LogRank", pmethod="Lau94")

temp$Group = as.factor(ifelse(temp$CE9>tempsplit$estimate,paste0(nazwa,">",round(tempsplit$estimate,2)),paste0(nazwa,"<=",round(tempsplit$estimate,2))))

s1 <- survfit(Surv(OS_TIME, OS) ~ Group, data = temp)
s1 %>% tbl_survfit(times = c(1,2,5), label_header = "**{time}-year survival**") %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0("Files/OStime_Survival_vs_",nazwa,".docx"))
png(paste0("Plots/OS_vs_",nazwa,".png"), width = 2000, height = 1500, res = 300)
ggsurvplot(s1, pval = T, risk.table = T, legend.labs = levels(as.factor(temp$Group)), break.time.by = 12, xlim=c(0,60))
dev.off()
m = coxph(Surv(OS_TIME, OS)  ~ Group, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0("Files/UnivariableOS_",nazwa,".docx"))

# Heatmap ----
library(ComplexHeatmap)



# Calculate TNM and platin adjusted survival ----
temp = merged_d

tbl_uv <-
  tbl_uvregression(
    dplyr::select(merged_d, -ID),
    method = coxph,
    y = Surv(OS_TIME, OS),
    exponentiate = TRUE,
    pvalue_fun = function(x) style_pvalue(x, digits = 2)
  )
tbl_uv %>% add_global_p() %>% add_q() %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0("AdjustedOS/Univariable_analysis.docx"))

colnames(merged_d)
m = coxph(Surv(OS_TIME, OS) ~ Cohort + Age_at_diagnosis + T + N + M + Platin, data = temp)
summary(m)

m2 = step(m)
m2

sink("AdjustedOS/Adjustment_model_stats.txt")
summary(m2)
sink()

predict(m2, type = "terms")
m2_survfit <- survfit(m2, newdata = temp)

merged_d$OS_TIME_ADJUSTED = summary(m2_survfit)$table[,'median']
merged_d$OS_TIME_ADJUSTED_RMEAN = summary(m2_survfit)$table[,'rmean']

tbl_uv <-
  tbl_uvregression(
    dplyr::select(merged_d, -ID),
    method = coxph,
    y = Surv(OS_TIME_ADJUSTED, OS),
    exponentiate = TRUE,
    pvalue_fun = function(x) style_pvalue(x, digits = 2)
  )
tbl_uv %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0("AdjustedOS/Univariable_analysis_adjusted.docx"))
tbl_uv

tbl_uv <-
  tbl_uvregression(
    dplyr::select(merged_d, -ID),
    method = coxph,
    y = Surv(OS_TIME_ADJUSTED_RMEAN, OS),
    exponentiate = TRUE,
    pvalue_fun = function(x) style_pvalue(x, digits = 2)
  )
tbl_uv %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0("AdjustedOS/Univariable_analysis_adjusted_rmean.docx"))
tbl_uv


# Intragroup associations ----
immuoparameters = colnames(merged_d)[c(immune_cols, ecotype_cols, consensus_cols)]


# Metaanalysis ----
setwd("H:/My Drive/PostCDDP/ANALYSIS/v4_032024/Metanalysis_results/General")
#immuoparameters = colnames(merged_d)[15:length(colnames(merged_d))-2]
#immuoparameters = immuoparameters[immuoparameters!="Operable"]
if(!dir.exists("reports")) { dir.create("reports") }
if(!dir.exists("forests")) { dir.create("forests") }
merged_d$Platin = as.factor(merged_d$Platin)
str(merged_d)
library(plyr)

kohorty = as.character(unique(merged_d$Cohort))
adjustedHRs = data.frame()

for(j in 1:length(kohorty)){
  
  kohorta = kohorty[j]
  temp = merged_d[merged_d$Cohort == kohorta,]
  #wyniki = list()
  wynikidf = data.frame()
  
  for(i in 1:length(immuoparameters)){
    try({
      selected = immuoparameters[i]
      tempxx = dplyr::select(temp, OS_TIME, OS, Age_at_diagnosis, T, M, Platin, all_of(selected))
      tempxx = tempxx[complete.cases(tempxx),]
      
      f = as.formula(paste0("Surv(OS_TIME, OS) ~ Age_at_diagnosis + 
    T + M + ", selected))
      m = coxph(f , data = tempxx)
      summary(m)
      
      wynikidf[i,"parameter"] = selected
      wynikidf[i,"loglik"] = logLik(m)
      wynikidf[i,"HR"] = summary(m)$conf.int[4,1]
      wynikidf[i,"LowerCI"] = summary(m)$conf.int[4,3]
      wynikidf[i,"UpperCI"] = summary(m)$conf.int[4,4]
      wynikidf[i,"cindex"] = summary(m)$concordance["C"]
      wynikidf[i,"pval"] = summary(m)$coefficients[4,5]
      wynikidf[i,"cohort"] = kohorta
      wynikidf[i,"n"] = sum(temp$Cohort == kohorta)
    })
    
    
    #a = tbl_regression(m, exponentiate = T)
    #wyniki = c(wyniki, list(a))
  }
  
  adjustedHRs = rbind.fill(adjustedHRs, wynikidf)
}
data.table::fwrite(adjustedHRs, "adjustedHRs.csv")

library(meta)
library(metafor)

wyniki_metaanaliz = data.frame()

for(i in 1:length(immuoparameters)) {
  parametr = immuoparameters[i]
  temp = dplyr::filter(adjustedHRs, parameter == parametr)
  temp$logHR = log(temp$HR)
  temp$logUB = log(temp$UpperCI)
  temp$logLB = log(temp$LowerCI)
  temp$selogHR =(temp$logUB - temp$logLB)/(2 * 1.96) 
  
  metaanalysis = NULL
  metaanalysis <- metagen(logHR, selogHR, studlab = cohort, data = temp, sm = "HR")
  
  if(!is.null(metaanalysis)) {
    sink(paste0("reports/",parametr,".txt"))
    print(summary(metaanalysis))
    sink()
    
    # summary(metaanalysis)
    # meta::forest(metaanalysis, JAMA.pval = T)
    
    wyniki_metaanaliz[i,"parameter"] = parametr
    wyniki_metaanaliz[i,"HR (fixed)"] =  exp(metaanalysis$TE.fixed)
    wyniki_metaanaliz[i,"Lower CI (fixed)"] =  exp(metaanalysis$lower.fixed)
    wyniki_metaanaliz[i,"Upper CI (fixed)"] =  exp(metaanalysis$upper.fixed)
    wyniki_metaanaliz[i,"p (fixed)"] =  metaanalysis$pval.fixed
    wyniki_metaanaliz[i,"I2"] =  metaanalysis$I2
    wyniki_metaanaliz[i,"HR (random)"] =  exp(metaanalysis$TE.random)
    wyniki_metaanaliz[i,"Lower CI (random)"] =  exp(metaanalysis$lower.random)
    wyniki_metaanaliz[i,"Upper CI (random)"] =  exp(metaanalysis$upper.random)
    wyniki_metaanaliz[i,"p (random)"] =  metaanalysis$pval.random
    
    pdf(paste0("forests/",parametr,".pdf"), width = 20)
    meta::forest(metaanalysis, title = parametr)
    dev.off()
    
    
    try({
      a = metabias(metaanalysis, method="rank", k.min=2)
      wyniki_metaanaliz[i,"BeggMazumdar p"] = a$pval
    })
    
    try({
      a = metabias(metaanalysis, method="linreg", k.min=2)
      wyniki_metaanaliz[i,"Egger p"] = a$pval
    })
    
  }
  
}

data.table::fwrite(wyniki_metaanaliz, "metaanalyses.csv")


# by Cisplatin ----
setwd("H:/My Drive/PostCDDP/ANALYSIS/v4_032024/Metanalysis_results/Cisplatin_vs_Carboplatin")
if(!dir.exists("reports")) { dir.create("reports") }
if(!dir.exists("forests")) { dir.create("forests") }

table(merged_d$Cohort)
table(merged_d$Platin)
table(merged_d$Cohort, merged_d$Platin)
#chisq.test(table(metadane$Cohort, dane$Cisplatin))

meta = paste0(merged_d$Cohort, " + Platin=",merged_d$Platin)
kohorty = unique(meta)

# Quick checks
temp = merged_d
# merged_d$meta = meta
# ggboxplot(temp, "meta", "B.cells.memory", add = "jitter") + stat_compare_means() + ylab("B cells memory") + xlab("")
# ggboxplot(temp, "meta", "B_cells", add = "jitter") + stat_compare_means() + ylab("B cells") + xlab("")
# ggboxplot(temp, "meta", "B_and_plasma_cells", add = "jitter") + stat_compare_means() + ylab("B cells (incl. plasma cells)") + xlab("")

adjustedHRs = data.frame()

for(j in 1:length(kohorty)){
  
  kohorta = kohorty[j]
  temp = merged_d[merged_d$meta == kohorta,]
  #wyniki = list()
  wynikidf = data.frame()
  
  for(i in 1:length(immuoparameters)){
    try({
      selected = immuoparameters[i]
      tempxx = dplyr::select(temp, OS_TIME, OS, Age_at_diagnosis, T, M, Platin, all_of(selected))
      tempxx = tempxx[complete.cases(tempxx),]
      
      f = as.formula(paste0("Surv(OS_TIME, OS) ~ Age_at_diagnosis + 
    T + M + ", selected))
      m = coxph(f , data = tempxx)
      summary(m)
      
      wynikidf[i,"parameter"] = selected
      wynikidf[i,"loglik"] = logLik(m)
      wynikidf[i,"HR"] = summary(m)$conf.int[4,1]
      wynikidf[i,"LowerCI"] = summary(m)$conf.int[4,3]
      wynikidf[i,"UpperCI"] = summary(m)$conf.int[4,4]
      wynikidf[i,"cindex"] = summary(m)$concordance["C"]
      wynikidf[i,"pval"] = summary(m)$coefficients[4,5]
      wynikidf[i,"cohort"] = kohorta
      wynikidf[i,"n"] = sum(temp$Cohort == kohorta)
    })
    
    
    #a = tbl_regression(m, exponentiate = T)
    #wyniki = c(wyniki, list(a))
  }
  
  adjustedHRs = rbind.fill(adjustedHRs, wynikidf)
}
data.table::fwrite(adjustedHRs, "adjustedHRs.csv")

library(meta)
library(metafor)

wyniki_metaanaliz = data.frame()
for(i in 1:length(immuoparameters)) {
  a = try({
    parametr = immuoparameters[i]
    temp = dplyr::filter(adjustedHRs, parameter == parametr)
    temp$logHR = log(temp$HR)
    temp$logUB = log(temp$UpperCI)
    temp$logLB = log(temp$LowerCI)
    temp$selogHR =(temp$logUB - temp$logLB)/(2 * 1.96) 
    
    temp = temp[complete.cases(temp),]
    # temp
    
    metaanalysis <- metagen(logHR, selogHR, studlab = cohort, data = temp, sm = "HR")
    
    
    # groups
    temp$subgroups = unlist(lapply(strsplit(temp$cohort, split = " \\+ "), function(x) x[2]))
    subgroups = unlist(lapply(strsplit(temp$cohort, split = " \\+ "), function(x) x[2]))
    metaregression = metareg(metaanalysis, ~subgroups)
    metaanalysis = update(metaanalysis, byvar = temp$subgroups)
    
    sink(paste0("reports/",parametr,".txt"))
    print(summary(metaanalysis))
    print(summary(metaregression))
    sink()
    
    # summary(metaanalysis)
    # meta::forest(metaanalysis, JAMA.pval = T)
    
    wyniki_metaanaliz[i,"parameter"] = parametr
    wyniki_metaanaliz[i,"HR (fixed)"] =  exp(metaanalysis$TE.fixed)
    wyniki_metaanaliz[i,"Lower CI (fixed)"] =  exp(metaanalysis$lower.fixed)
    wyniki_metaanaliz[i,"Upper CI (fixed)"] =  exp(metaanalysis$upper.fixed)
    wyniki_metaanaliz[i,"p (fixed)"] =  metaanalysis$pval.fixed
    wyniki_metaanaliz[i,"I2"] =  metaanalysis$I2
    wyniki_metaanaliz[i,"HR (random)"] =  exp(metaanalysis$TE.random)
    wyniki_metaanaliz[i,"Lower CI (random)"] =  exp(metaanalysis$lower.random)
    wyniki_metaanaliz[i,"Upper CI (random)"] =  exp(metaanalysis$upper.random)
    wyniki_metaanaliz[i,"p (random)"] =  metaanalysis$pval.random
    
    wyniki_metaanaliz[i,"Common between p"] =  metaanalysis$pval.Q.b.common
    wyniki_metaanaliz[i,"Random between p"] =  metaanalysis$pval.Q.b.random
    
    wyniki_metaanaliz[i,"Metaregression QM"] =  metaregression$QM
    wyniki_metaanaliz[i,"Metaregression p"] =  metaregression$QMp
    
    pdf(paste0("forests/",parametr,".pdf"), width = 20, height = 20)
    meta::forest(metaanalysis)
    #meta::forest(metaregression)
    bubble(metaregression, studlab = T)
    
    dev.off()
    
    # a =a$pval metabias(metaanalysis, method="rank", k.min=6)
    # wyniki_metaanaliz[i,"BeggMazumdar p"] = a$pval
    # 
    # a = metabias(metaanalysis, method="linreg", k.min=6)
    # wyniki_metaanaliz[i,"Egger p"] = 
  })
  if(class(a) == "try-error"){
    file.remove(paste0("forests/",parametr,".pdf"))
  }
}
data.table::fwrite(wyniki_metaanaliz, "metaanalyses.csv")

# B-cell analysis ----
#setwd("/drive/My Drive/PostCDDP/ANALYSIS/v4_032024/B-cell_analysis")

library(maxstat)
temp = merged_d
nazwa = "B-cells memory"
temp$B.cells.memory
tempsplit = maxstat.test(Surv(OS_TIME, OS) ~ B.cells.memory,  data=temp, smethod="LogRank", pmethod="Lau94")

temp$Group = as.factor(ifelse(temp$B.cells.memory>tempsplit$estimate,paste0(nazwa,">",round(tempsplit$estimate,2)),paste0(nazwa,"<=",round(tempsplit$estimate,2))))
table(temp$Group)
merged_d$B.cells.memory_groups = temp$Group

s1 <- survfit(Surv(OS_TIME, OS) ~ Group, data = temp)
#s1 %>% tbl_survfit(times = c(12,24,60,120), label_header = "**{time}-year survival**") %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa, " unadjusted high_vs_low - Survival Table.docx"))
png(paste0(nazwa," unadjusted high_vs_low - KM.png"), width = 2000, height = 1500, res = 300)
ggsurvplot(s1, pval = T, risk.table = T, legend.labs = levels(as.factor(temp$Group)), break.time.by = 12, xlim=c(0,120))
dev.off()
m = coxph(Surv(OS_TIME, OS) ~ Group, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa," unadjusted high_vs_low - HR.docx"))

##
library(maxstat)
temp = merged_d
nazwa = "CE10"
temp$CE10
tempsplit = maxstat.test(Surv(OS_TIME, OS) ~CE10,  data=temp, smethod="LogRank", pmethod="Lau94")

temp$Group = as.factor(ifelse(temp$CE10>tempsplit$estimate,paste0(nazwa,">",round(tempsplit$estimate,2)),paste0(nazwa,"<=",round(tempsplit$estimate,2))))
table(temp$Group)
merged_d$CE10_groups = temp$Group

s1 <- survfit(Surv(OS_TIME, OS) ~ Group, data = temp)
#s1 %>% tbl_survfit(times = c(12,24,60,120), label_header = "**{time}-year survival**") %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa, " unadjusted high_vs_low - Survival Table.docx"))
png(paste0(nazwa," unadjusted high_vs_low - KM.png"), width = 2000, height = 1500, res = 300)
ggsurvplot(s1, pval = T, risk.table = T, legend.labs = levels(as.factor(temp$Group)), break.time.by = 12, xlim=c(0,120))
dev.off()
m = coxph(Surv(OS_TIME, OS) ~ Group, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa," unadjusted high_vs_low - HR.docx"))

##
library(maxstat)
temp = merged_d
nazwa = "CE2"
temp$CE2
tempsplit = maxstat.test(Surv(OS_TIME, OS) ~CE2,  data=temp, smethod="LogRank", pmethod="Lau94")

temp$Group = as.factor(ifelse(temp$CE2>tempsplit$estimate,paste0(nazwa,">",round(tempsplit$estimate,2)),paste0(nazwa,"<=",round(tempsplit$estimate,2))))
table(temp$Group)
merged_d$CE2_groups = temp$Group

s1 <- survfit(Surv(OS_TIME, OS) ~ Group, data = temp)
#s1 %>% tbl_survfit(times = c(12,24,60,120), label_header = "**{time}-year survival**") %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa, " unadjusted high_vs_low - Survival Table.docx"))
png(paste0(nazwa," unadjusted high_vs_low - KM.png"), width = 2000, height = 1500, res = 300)
ggsurvplot(s1, pval = T, risk.table = T, legend.labs = levels(as.factor(temp$Group)), break.time.by = 12, xlim=c(0,120))
dev.off()
m = coxph(Surv(OS_TIME, OS) ~ Group, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa," unadjusted high_vs_low - HR.docx"))

##
table(merged_d$B.cells.memory_groups, merged_d$Cohort)
chisq.test(table(merged_d$B.cells.memory_groups, merged_d$Cohort))

table(merged_d$B.cells.memory_groups, merged_d$CE10_groups)
chisq.test(table(merged_d$B.cells.memory_groups, merged_d$CE10_groups))
table(merged_d$B.cells.memory_groups, merged_d$CE2_groups)
chisq.test(table(merged_d$B.cells.memory_groups, merged_d$CE2_groups))

table(merged_d$B.cells.memory_groups, merged_d$Cisplatin)
chisq.test(table(merged_d$B.cells.memory_groups, merged_d$Cohort))

chisq.test(table(merged_d$B.cells.memory_groups, merged_d$Platin))
chisq.test(table(merged_d$B.cells.memory_groups, merged_d$T))
chisq.test(table(merged_d$B.cells.memory_groups, merged_d$N))
chisq.test(table(merged_d$B.cells.memory_groups, merged_d$M))
chisq.test(table(merged_d$B.cells.memory_groups, merged_d$`Carcinoma Ecotype`))

table(merged_d$B.cells.memory_groups, merged_d$consensusClass)
chisq.test(table(merged_d$B.cells.memory_groups, merged_d$consensusClass))

chisq.test(table(merged_d$B.cells.memory_groups, merged_d$Baylor.subtype))

merged_d$Myeloid
a = ggboxplot(merged_d, "B.cells.memory_groups", "Age_at_diagnosis", add = "jitter") + stat_compare_means(method = "t.test") + xlab("B-cell memory groups") + ylab("Age at diagnosis")
b = ggboxplot(merged_d, "B.cells.memory_groups", "T.cells.CD8", add = "jitter") + stat_compare_means(method = "t.test") + xlab("B-cell memory groups") + ylab("T-cells CD8")
c = ggboxplot(merged_d, "B.cells.memory_groups", "Tcells_CD4", add = "jitter") + stat_compare_means(method = "t.test") + xlab("B-cell memory groups") + ylab("T-cells CD4")
d = ggboxplot(merged_d, "B.cells.memory_groups", "Monocytes_derivatives", add = "jitter") + stat_compare_means(method = "t.test") + xlab("B-cell memory groups") + ylab("Monocytes derivatives")
e = ggboxplot(merged_d, "B.cells.memory_groups", "Myeloid", add = "jitter") + stat_compare_means(method = "t.test") + xlab("B-cell memory groups") + ylab("Myeloid lineage")
f = ggboxplot(merged_d, "B.cells.memory_groups", "Lymphoid", add = "jitter") + stat_compare_means(method = "t.test") + xlab("B-cell memory groups") + ylab("Lymphoid lineage")
ggarrange(plotlist = list(a,b,c,d,e,f), nrow = 3, ncol =2)
####
temp$B.cells.memory_groups
s1 <- survfit(Surv(OS_TIME, OS) ~ CE10_groups + Platin, data = temp)
e = ggsurvplot(s1, pval = T, risk.table = T, break.time.by=12, xlim = c(0,84))
e

ggsave(file = "Plots/OS_vs_B.cells.memory_groups + Platin.png", plot= e, width = 8, height = 6, dpi = 300)




###


tempsplit = maxstat.test(Surv(OS_TIME_ADJUSTED_RMEAN, OS) ~ B.cells.memory,  data=temp, smethod="LogRank", pmethod="Lau94")

temp$Group = as.factor(ifelse(temp$B.cells.memory>tempsplit$estimate,paste0(nazwa,">",round(tempsplit$estimate,2)),paste0(nazwa,"<=",round(tempsplit$estimate,2))))
table(temp$Group)
merged_d$B.cells.memory_groups_radjusted = temp$Group
temp$Group_adjusted = temp$Group
s1 <- survfit(Surv(OS_TIME_ADJUSTED_RMEAN, OS) ~ Group, data = temp)
s1 %>% tbl_survfit(times = c(12,24,60,120), label_header = "**{time}-year survival**") %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa, " rmean-adjusted high_vs_low - Survival Table.docx"))
png(paste0(nazwa," rmean-adjusted high_vs_low - KM.png"), width = 2000, height = 1500, res = 300)
ggsurvplot(s1, pval = T, risk.table = T, legend.labs = levels(as.factor(temp$Group)), break.time.by = 12, xlim=c(0,120))
dev.off()
m = coxph(Surv(OS_TIME_ADJUSTED_RMEAN, OS) ~ Group, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa," rmean-adjusted high_vs_low - HR.docx"))


library(maxstat)
temp = merged_d
nazwa = "B-cells"
temp$B_cells
tempsplit = maxstat.test(Surv(OS_TIME, OS) ~ B_cells,  data=temp, smethod="LogRank", pmethod="Lau94")

temp$Group = as.factor(ifelse(temp$B_cells>tempsplit$estimate,paste0(nazwa,">",round(tempsplit$estimate,2)),paste0(nazwa,"<=",round(tempsplit$estimate,2))))
table(temp$Group)
merged_d$B_cells_groups = temp$Group
temp$Group_unadjusted = temp$Group

s1 <- survfit(Surv(OS_TIME, OS) ~ Group, data = temp)
s1 %>% tbl_survfit(times = c(12,24,60,120), label_header = "**{time}-year survival**") %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa, " unadjusted high_vs_low - Survival Table.docx"))
png(paste0(nazwa," unadjusted high_vs_low - KM.png"), width = 2000, height = 1500, res = 300)
ggsurvplot(s1, pval = T, risk.table = T, legend.labs = levels(as.factor(temp$Group)), break.time.by = 12, xlim=c(0,120))
dev.off()
m = coxph(Surv(OS_TIME, OS) ~ Group, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa," unadjusted high_vs_low - HR.docx"))

tempsplit = maxstat.test(Surv(OS_TIME_ADJUSTED, OS) ~ B_cells,  data=temp, smethod="LogRank", pmethod="Lau94")

temp$Group = as.factor(ifelse(temp$B.cells.memory>tempsplit$estimate,paste0(nazwa,">",round(tempsplit$estimate,2)),paste0(nazwa,"<=",round(tempsplit$estimate,2))))
table(temp$Group)
merged_d$B_cells_groups_adjusted = temp$Group
s1 <- survfit(Surv(OS_TIME_ADJUSTED, OS) ~ Group, data = temp)
s1 %>% tbl_survfit(times = c(12,24,60,120), label_header = "**{time}-year survival**") %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa, " adjusted high_vs_low - Survival Table.docx"))
png(paste0(nazwa," adjusted high_vs_low - KM.png"), width = 2000, height = 1500, res = 300)
ggsurvplot(s1, pval = T, risk.table = T, legend.labs = levels(as.factor(temp$Group)), break.time.by = 12, xlim=c(0,120))
dev.off()
m = coxph(Surv(OS_TIME_ADJUSTED, OS) ~ Group, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa," adjusted high_vs_low - HR.docx"))

tempsplit = maxstat.test(Surv(OS_TIME_ADJUSTED_RMEAN, OS) ~ B_cells,  data=temp, smethod="LogRank", pmethod="Lau94")

temp$Group = as.factor(ifelse(temp$B_cells>tempsplit$estimate,paste0(nazwa,">",round(tempsplit$estimate,2)),paste0(nazwa,"<=",round(tempsplit$estimate,2))))
table(temp$Group)
merged_d$B_cells_groups_radjusted = temp$Group
temp$Group_adjusted = temp$Group
s1 <- survfit(Surv(OS_TIME_ADJUSTED_RMEAN, OS) ~ Group, data = temp)
s1 %>% tbl_survfit(times = c(12,24,60,120), label_header = "**{time}-year survival**") %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa, " rmean-adjusted high_vs_low - Survival Table.docx"))
png(paste0(nazwa," rmean-adjusted high_vs_low - KM.png"), width = 2000, height = 1500, res = 300)
ggsurvplot(s1, pval = T, risk.table = T, legend.labs = levels(as.factor(temp$Group)), break.time.by = 12, xlim=c(0,120))
dev.off()
m = coxph(Surv(OS_TIME_ADJUSTED_RMEAN, OS) ~ Group, data = temp)
summary(m)
tbl_regression(m, exponentiate = T) %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa," rmean-adjusted high_vs_low - HR.docx"))



## Plot of individuals ----
temp$OS = as.factor(temp$OS)
table(temp$OS)
temp$Survival = ifelse(temp$OS == 1, "Death", "Survival")
a = temp %>%
  ggplot(aes(x = ID, y = OS_TIME)) +
  geom_linerange(aes(ymin = 0, ymax = OS_TIME, color = Cohort) ) +
  geom_point(aes(shape = Survival), stroke = 1, cex = 1) +
  scale_shape_manual(values = c(0, NA)) +
  labs(y = "Time (months)", x = "Subjects") + coord_flip() + theme_classic() + theme(axis.text.y = element_blank())
a
ggsave("Unadjusted profiles.png", plot = a, dpi = 300, width = 5, height = 10)


temp$OS = as.factor(temp$OS)
a = temp %>%
  ggplot(aes(x = ID, y = OS_TIME_ADJUSTED_RMEAN)) +
  geom_linerange(aes(ymin = 0, ymax = OS_TIME_ADJUSTED_RMEAN, color = Cohort) ) +
  geom_point(aes(shape = OS), stroke = 1, cex = 1) +
  scale_shape_manual(values = c(NA, 0)) +
  labs(y = "Time (months)", x = "Subjects") + coord_flip() + theme_classic() + theme(axis.text.y = element_blank())
ggsave("Rmean-adjusted profiles.png", plot = a, dpi = 300, width = 5, height = 10)

# Selected signatures ----
setwd("H:/My Drive/PostCDDP/ANALYSIS/v4_032024")
library(dplyr)
TMM_HMGREEK = data.table::fread("../v1_Dec2022/HM+GREEK/Data/TMM.csv") %>% as.data.frame()
rownames(TMM_HMGREEK) = TMM_HMGREEK$Gene
TMM_HMGREEK = t(as.matrix(TMM_HMGREEK))
colnames(TMM_HMGREEK)
rownames(TMM_HMGREEK)
TMM_HMGREEK = TMM_HMGREEK[-1,]
TMM_HMGREEK = as.data.frame(TMM_HMGREEK)
TMM_HMGREEK$ID_for_matching = rownames(TMM_HMGREEK)

TMM_TABER = data.table::fread("../v1_Dec2022/TABER/Data/TMM.csv") %>% as.data.frame()
rownames(TMM_TABER) = TMM_TABER$Gene
TMM_TABER = t(as.matrix(TMM_TABER))
colnames(TMM_TABER)
rownames(TMM_TABER)
TMM_TABER = TMM_TABER[-1,]
TMM_TABER = as.data.frame(TMM_TABER)
TMM_TABER$ID_for_matching = rownames(TMM_TABER)

library(readr)
merged_d <- read_csv("Files/merged_d.csv")
merged_d$ID_for_matching = make.names(merged_d$ID)

library(plyr); library(dplyr)
TMM = rbind.fill(TMM_HMGREEK, TMM_TABER)

temp = dplyr::left_join(merged_d, TMM, by = "ID_for_matching", suffix = c("","_v2"))
all = as.data.frame(temp)

setwd("SelectedSignatures")
data.table::fwrite(temp, "metadane+TMM.csv.gz")

library(maxstat)
library(flextable)
library(gtsummary)
library(ggpubr)
library(survminer)

## B cell-related gene (BCR) signature consisting of nine cytokine signalling genes (Br J Cancer 2022) ----
nazwa = "Br J Cancer 2022"
signature = c("CXCR6", "IL18RAP", "LCK", "IL2RG", "CXCL13", "PSMB10", "TNFRSF14", "BATF", "TNFRSF4")
x = dplyr::select(all, all_of(signature)) %>% as.matrix()
mode(x) = "numeric"
all$bcell_signature1 = rowMeans(x)

cor.test(all$bcell_signature1, all$B.cells.memory)
#source("https://raw.githubusercontent.com/kstawiski/OmicSelector/master/R/OmicSelector_correlation_plot.R")

png(paste0(nazwa," high_vs_low - correlation.png"), width = 2000, height = 1500, res = 300)
OmicSelector_correlation_plot(all$bcell_signature1, all$B.cells.memory, labvar1 = "B cell-related nine cytokine signalling genes (Br J Cancer 2022)", labvar2 = "B-cell memory (CIBERSORTx)", yx = F, title ="", metoda = "spearman")
dev.off()

tempsplit = maxstat.test(Surv(OS_TIME, OS) ~ bcell_signature1,  data=all, smethod="LogRank", pmethod="Lau94")
all$Group = as.factor(ifelse(all$bcell_signature1>tempsplit$estimate,"High score","Low score"))
table(all$Group)
all$bcell_signature1_group = all$Group
s1 <- survfit(Surv(OS_TIME, OS) ~ Group, data = as.data.frame(all))
s1 %>% tbl_survfit(times = c(12,24,60,120), label_header = "**{time}-month survival**") %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa, " high_vs_low - Survival Table.docx"))
png(paste0(nazwa," high_vs_low - KM.png"), width = 2000, height = 1500, res = 300)
ggsurvplot(s1, pval = T, risk.table = T, legend.labs = levels(as.factor(all$Group)), break.time.by = 12, xlim=c(0,120), data = all, title = nazwa) 
dev.off()
m = coxph(Surv(OS_TIME, OS) ~ Group, data = all)
summary(m)
tbl_regression(m, exponentiate = T) %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa," high_vs_low - HR.docx"))

## B lineage-associated risk signature predicting prognosis of patients with lung adenocarcinoma (BMC Cancer 2021) ----
# -0.16936*FCER2;-0.35238*CR2;-0.16979*FCRLA;-0.27845*BLK;-0.32492*MS4A1;-0.21464*FCRL1;-0.01268*ALB;-0.03807*KRT20;-0.63325*CD19;-0.17511*FDCSP;-0.08366*CNR2;-0.05848*GH1;-0.47873*TNFRSF13B
nazwa = "BMC Cancer 2021"
all$GH1[is.na(all$GH1)] = -2.57271994
all$bcell_signature2 = -0.16936*as.numeric(all$FCER2) + -0.35238*as.numeric(all$CR2) + -0.16979*as.numeric(all$FCRL1) + -0.27845*as.numeric(all$BLK) + -0.32492*as.numeric(all$MS4A1) + -0.21464*as.numeric(all$FCRL1) + -0.01268*as.numeric(all$ALB) + -0.03807*as.numeric(all$KRT20) + -0.63325*as.numeric(all$CD19) + -0.17511*as.numeric(all$FDCSP) + -0.08366*as.numeric(all$CNR2) + -0.05848*as.numeric(all$GH1) + -0.47873*as.numeric(all$TNFRSF13B)


cor.test(all$bcell_signature2, all$B.cells.memory)
#source("https://raw.githubusercontent.com/kstawiski/OmicSelector/master/R/OmicSelector_correlation_plot.R")

png(paste0(nazwa," high_vs_low - correlation.png"), width = 2000, height = 1500, res = 300)
OmicSelector_correlation_plot(-all$bcell_signature2, all$B.cells.memory, labvar1 = "Inv. B lineage-associated risk signature, lung adenocarcinoma (BMC Cancer 2021)", labvar2 = "B-cell memory (CIBERSORTx)", yx = F, title ="", metoda = "spearman")
dev.off()

tempsplit = maxstat.test(Surv(OS_TIME, OS) ~ bcell_signature2,  data=all, smethod="LogRank", pmethod="Lau94")
all$Group = as.factor(ifelse(all$bcell_signature2>tempsplit$estimate,"High score","Low score"))
table(all$Group)
all$bcell_signature2_group = all$Group
s1 <- survfit(Surv(OS_TIME, OS) ~ Group, data = as.data.frame(all))
s1 %>% tbl_survfit(times = c(12,24,60,120), label_header = "**{time}-month survival**") %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa, " high_vs_low - Survival Table.docx"))
png(paste0(nazwa," high_vs_low - KM.png"), width = 2000, height = 1500, res = 300)
ggsurvplot(s1, pval = T, risk.table = T, legend.labs = levels(as.factor(all$Group)), break.time.by = 12, xlim=c(0,120), data = all, title = nazwa) 
dev.off()
m = coxph(Surv(OS_TIME, OS) ~ Group, data = all)
summary(m)
tbl_regression(m, exponentiate = T) %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa," high_vs_low - HR.docx"))



## Identification of a Five-Gene Prognostic Signature Related to B Cells Infiltration in Pancreatic Adenocarcinoma  Int J Gen Med. 2021 Aug 30;14:5051-5068. doi: 10.2147/IJGM.S324432. eCollection 2021. ----
# Risk score = (−0.502) * ARID5A expression level + 1.226 * CLEC2B expression level + (−0.419) * MICAL1 expression level + (−0.309) * MZB1 expression level + (−0.757) * RAPGEF1 expression level
nazwa = "Int J Gen Med. 2021"
all$bcell_signature3 = -0.502*as.numeric(all$ARID5A) + 1.226*as.numeric(all$CLEC2B) + -0.419*as.numeric(all$MICAL1) + -0.309*as.numeric(all$MZB1) + -0.757*as.numeric(all$RAPGEF1)

cor.test(-all$bcell_signature3, all$B.cells.memory)
all$bcell_signature3 = -all$bcell_signature3

source("https://raw.githubusercontent.com/kstawiski/OmicSelector/master/R/OmicSelector_correlation_plot.R")

png(paste0(nazwa," high_vs_low - correlation.png"), width = 2000, height = 1500, res = 300)
OmicSelector_correlation_plot(all$bcell_signature3, all$B.cells.memory, labvar1 = "B Cells Infiltration in Pancreatic Adenocarcinoma (Int J Gen Med. 2021)", labvar2 = "B-cell memory (CIBERSORTx)", yx = F, title ="", metoda = "spearman")
dev.off()

tempsplit = maxstat.test(Surv(OS_TIME, OS) ~ bcell_signature3,  data=all, smethod="LogRank", pmethod="Lau94")
all$Group = as.factor(ifelse(all$bcell_signature3>tempsplit$estimate,"High score","Low score"))
table(all$Group)
all$bcell_signature3_group = all$Group
s1 <- survfit(Surv(OS_TIME, OS) ~ Group, data = as.data.frame(all))
s1 %>% tbl_survfit(times = c(12,24,60,120), label_header = "**{time}-month survival**") %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa, " high_vs_low - Survival Table.docx"))
png(paste0(nazwa," high_vs_low - KM.png"), width = 2000, height = 1500, res = 300)
ggsurvplot(s1, pval = T, risk.table = T, legend.labs = levels(as.factor(all$Group)), break.time.by = 12, xlim=c(0,120), data = all, title = nazwa) 
dev.off()
m = coxph(Surv(OS_TIME, OS) ~ Group, data = all)
summary(m)
tbl_regression(m, exponentiate = T) %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa," high_vs_low - HR.docx"))


## Gene signature based on B cell predicts clinical outcome of radiotherapy and immunotherapy for patients with lung adenocarcinoma (Cancer Med 2020) ----
#  -3.503*PARP15;-1.624*FADS3;-2.0592*RUBCNL;2.26*BIRC3;3.385*SP110;2.286*TLE1
nazwa = "Cancer Med 2020"
all$RUBCNL[is.na(all$RUBCNL)] = -2.57
all$bcell_signature4 = -3.503*as.numeric(all$PARP15) + -1.624*as.numeric(all$FADS3) + -2.0592*as.numeric(all$RUBCNL) + 2.26*as.numeric(all$BIRC3) + 3.385*as.numeric(all$SP110) + 2.286*as.numeric(all$TLE1)

cor.test(all$bcell_signature4, all$B.cells.memory)

source("https://raw.githubusercontent.com/kstawiski/OmicSelector/master/R/OmicSelector_correlation_plot.R")

png(paste0(nazwa," high_vs_low - correlation.png"), width = 2000, height = 1500, res = 300)
OmicSelector_correlation_plot(all$bcell_signature4, all$B.cells.memory, labvar1 = "Gene signature based on B cell, lung adenocarcinoma (Cancer Med 2020)", labvar2 = "B-cell memory (CIBERSORTx)", yx = F, title ="", metoda = "spearman")
dev.off()

tempsplit = maxstat.test(Surv(OS_TIME, OS) ~ bcell_signature4,  data=all, smethod="LogRank", pmethod="Lau94")
all$Group = as.factor(ifelse(all$bcell_signature4>tempsplit$estimate,"High score","Low score"))
table(all$Group)
all$bcell_signature4_group = all$Group
s1 <- survfit(Surv(OS_TIME, OS) ~ Group, data = as.data.frame(all))
s1 %>% tbl_survfit(times = c(12,24,60,120), label_header = "**{time}-month survival**") %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa, " high_vs_low - Survival Table.docx"))
png(paste0(nazwa," high_vs_low - KM.png"), width = 2000, height = 1500, res = 300)
ggsurvplot(s1, pval = T, risk.table = T, legend.labs = levels(as.factor(all$Group)), break.time.by = 12, xlim=c(0,120), data = all, title = nazwa) 
dev.off()
m = coxph(Surv(OS_TIME, OS) ~ Group, data = all)
summary(m)
tbl_regression(m, exponentiate = T) %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa," high_vs_low - HR.docx"))

# A 14-gene B-cell immune signature in early-stage triple-negative breast cancer (TNBC): a pooled analysis of seven studies. EBioMedicine 2024 Apr:102:105043 ----
# The IGG signature is composed of 14 genes implicated in maturation of T and B lymphocytes progenitors (IL2RG), CD4+ and B lymphocytes activation and survival (CD27, TNFRSF17, PIM2), B lymphocytes differentiation in germinal centers (POU2AF1), immunoglobulin production (CD79a, JCHAIN, IGKC, IGL, IGLV3-25), chemotaxis (CXCL8, NTN3), and regulation of B, T and NK lymphocytes activity (LAX1, HLA-C). In datasets where the expression of some IGG genes was missing, the signature score was calculated from the available genes (Figure S1). For consistency with previous studies, the signature was calculated as mean expression of genes.
nazwa = "EBioMedicine 2024"
signature = c("IL2RG", "CD27", "TNFRSF17", "PIM2", "POU2AF1", "CD79A", "JCHAIN", "IGKC", "CXCL8", "NTN3", "LAX1", "HLA-C")

x = dplyr::select(all, all_of(signature)) %>% as.matrix()
#naniar::vis_miss(as.data.frame(x))
#x = x[,colSums(is.na(x))==0]
mode(x) = "numeric"
all$bcell_signature5 = rowMeans(x, na.rm = T)

source("https://raw.githubusercontent.com/kstawiski/OmicSelector/master/R/OmicSelector_correlation_plot.R")

png(paste0(nazwa," high_vs_low - correlation.png"), width = 2000, height = 1500, res = 300)
OmicSelector_correlation_plot(all$bcell_signature5, all$B.cells.memory, labvar1 = "14-gene B-cell immune signature in early-stage TNBC (EBioMedicine 2024)", labvar2 = "B-cell memory (CIBERSORTx)", yx = F, title ="", metoda = "spearman")
dev.off()

tempsplit = maxstat.test(Surv(OS_TIME, OS) ~ bcell_signature5,  data=all, smethod="LogRank", pmethod="Lau94")
all$Group = as.factor(ifelse(all$bcell_signature5>=tempsplit$estimate,"High score","Low score"))
table(all$Group)
all$bcell_signature5_group = all$Group
s1 <- survfit(Surv(OS_TIME, OS) ~ Group, data = as.data.frame(all))
s1 %>% tbl_survfit(times = c(12,24,60,120), label_header = "**{time}-month survival**") %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa, " high_vs_low - Survival Table.docx"))
png(paste0(nazwa," high_vs_low - KM.png"), width = 2000, height = 1500, res = 300)
ggsurvplot(s1, pval = T, risk.table = T, legend.labs = levels(as.factor(all$Group)), break.time.by = 12, xlim=c(0,120), data = all, title = nazwa) 
dev.off()
m = coxph(Surv(OS_TIME, OS) ~ Group, data = all)
summary(m)
tbl_regression(m, exponentiate = T)
tbl_regression(m, exponentiate = T) %>% as_flex_table() %>% flextable::save_as_docx(., path=paste0(nazwa," high_vs_low - HR.docx"))

data.table::fwrite(all, "all.csv.gz")

## Assoc ----
all$bcell
all$B.cells.memory
tempsplit = maxstat.test(Surv(OS_TIME, OS) ~ B.cells.memory,  data=all, smethod="LogRank", pmethod="Lau94")

nazwa = "B-cells memory"
all$B_cell_memory_groups = as.factor(ifelse(all$B.cells.memory>tempsplit$estimate,paste0(nazwa,">",round(tempsplit$estimate,2)),paste0(nazwa,"<=",round(tempsplit$estimate,2))))
table(all$B_cell_memory_groups)

table(all$B_cell_memory_groups, all$bcell_signature1_group)
chisq.test(table(all$B_cell_memory_groups, all$bcell_signature1_group))
tbl_regression(coxph(Surv(OS_TIME, OS) ~ bcell_signature1, data = all), exponentiate = T)

table(all$B_cell_memory_groups, all$bcell_signature2_group)
chisq.test(table(all$B_cell_memory_groups, all$bcell_signature2_group))
tbl_regression(coxph(Surv(OS_TIME, OS) ~ bcell_signature2, data = all), exponentiate = T)

table(all$B_cell_memory_groups, all$bcell_signature3_group)
chisq.test(table(all$B_cell_memory_groups, all$bcell_signature3_group))
tbl_regression(coxph(Surv(OS_TIME, OS) ~ bcell_signature3, data = all), exponentiate = T)

table(all$B_cell_memory_groups, all$bcell_signature4_group)
chisq.test(table(all$B_cell_memory_groups, all$bcell_signature4_group))
tbl_regression(coxph(Surv(OS_TIME, OS) ~ bcell_signature4, data = all), exponentiate = T)

table(all$B_cell_memory_groups, all$bcell_signature5_group)
chisq.test(table(all$B_cell_memory_groups, all$bcell_signature5_group))
tbl_regression(coxph(Surv(OS_TIME, OS) ~ bcell_signature5, data = all), exponentiate = T)

# Meta-analysis of signatures ----
all$Platin

immuoparameters = colnames(all)[startsWith(colnames(all),"bcell_sig")]
if(!dir.exists("reports")) { dir.create("reports") }
if(!dir.exists("forests")) { dir.create("forests") }
all$Platin = as.factor(all$Platin)
library(plyr)

kohorty = as.character(unique(all$Cohort))
adjustedHRs = data.frame()

for(j in 1:length(kohorty)){
  
  kohorta = kohorty[j]
  temp = all[all$Cohort == kohorta,]
  #wyniki = list()
  wynikidf = data.frame()
  
  for(i in 1:length(immuoparameters)){
    try({
      selected = immuoparameters[i]
      tempxx = dplyr::select(temp, OS_TIME, OS, Age_at_diagnosis, T, M, Platin, all_of(selected))
      tempxx = tempxx[complete.cases(tempxx),]
      
      f = as.formula(paste0("Surv(OS_TIME, OS) ~ Age_at_diagnosis + 
    T + M + ", selected))
      m = coxph(f , data = tempxx)
      summary(m)
      
      wynikidf[i,"parameter"] = selected
      wynikidf[i,"loglik"] = logLik(m)
      wynikidf[i,"HR"] = summary(m)$conf.int[4,1]
      wynikidf[i,"LowerCI"] = summary(m)$conf.int[4,3]
      wynikidf[i,"UpperCI"] = summary(m)$conf.int[4,4]
      wynikidf[i,"cindex"] = summary(m)$concordance["C"]
      wynikidf[i,"pval"] = summary(m)$coefficients[4,5]
      wynikidf[i,"cohort"] = kohorta
      wynikidf[i,"n"] = sum(temp$Cohort == kohorta)
    })
    
    
    #a = tbl_regression(m, exponentiate = T)
    #wyniki = c(wyniki, list(a))
  }
  
  adjustedHRs = rbind.fill(adjustedHRs, wynikidf)
}
data.table::fwrite(adjustedHRs, "adjustedHRs.csv")

library(meta)
library(metafor)

wyniki_metaanaliz = data.frame()

for(i in 1:length(immuoparameters)) {
  parametr = immuoparameters[i]
  temp = dplyr::filter(adjustedHRs, parameter == parametr)
  temp$logHR = log(temp$HR)
  temp$logUB = log(temp$UpperCI)
  temp$logLB = log(temp$LowerCI)
  temp$selogHR =(temp$logUB - temp$logLB)/(2 * 1.96) 
  
  metaanalysis = NULL
  metaanalysis <- metagen(logHR, selogHR, studlab = cohort, data = temp, sm = "HR")
  
  if(!is.null(metaanalysis)) {
    sink(paste0("reports/",parametr,".txt"))
    print(summary(metaanalysis))
    sink()
    
    # summary(metaanalysis)
    # meta::forest(metaanalysis, JAMA.pval = T)
    
    wyniki_metaanaliz[i,"parameter"] = parametr
    wyniki_metaanaliz[i,"HR (fixed)"] =  exp(metaanalysis$TE.fixed)
    wyniki_metaanaliz[i,"Lower CI (fixed)"] =  exp(metaanalysis$lower.fixed)
    wyniki_metaanaliz[i,"Upper CI (fixed)"] =  exp(metaanalysis$upper.fixed)
    wyniki_metaanaliz[i,"p (fixed)"] =  metaanalysis$pval.fixed
    wyniki_metaanaliz[i,"I2"] =  metaanalysis$I2
    wyniki_metaanaliz[i,"HR (random)"] =  exp(metaanalysis$TE.random)
    wyniki_metaanaliz[i,"Lower CI (random)"] =  exp(metaanalysis$lower.random)
    wyniki_metaanaliz[i,"Upper CI (random)"] =  exp(metaanalysis$upper.random)
    wyniki_metaanaliz[i,"p (random)"] =  metaanalysis$pval.random
    
    pdf(paste0("forests/",parametr,".pdf"), width = 20)
    meta::forest(metaanalysis, title = parametr)
    dev.off()
    
    
    try({
      a = metabias(metaanalysis, method="rank", k.min=2)
      wyniki_metaanaliz[i,"BeggMazumdar p"] = a$pval
    })
    
    try({
      a = metabias(metaanalysis, method="linreg", k.min=2)
      wyniki_metaanaliz[i,"Egger p"] = a$pval
    })
    
  }
  
}

wyniki_metaanaliz
data.table::fwrite(wyniki_metaanaliz, "metaanalyses.csv")


if(!dir.exists("byCisplatin")) { dir.create("byCisplatin") }
setwd("/drive/My Drive/PostCDDP/ANALYSIS/v4_032024/SelectedSignatures/byCisplatin")

if(!dir.exists("reports")) { dir.create("reports") }
if(!dir.exists("forests")) { dir.create("forests") }

table(all$Cohort)
table(all$Platin)
table(all$Cohort, all$Platin)
#chisq.test(table(metadane$Cohort, dane$Cisplatin))

meta = paste0(all$Cohort, " + Platin=",all$Platin)
kohorty = unique(meta)

# Quick checks
temp = all
all$meta = meta
all = left_join(all, merged_d, by = "ID", suffix = c("","v2"))
all$meta
ggboxplot(temp, "meta", "B.cells.memory", add = "jitter") + stat_compare_means() + ylab("B cells memory") + xlab("")
ggboxplot(temp, "meta", "B_cells", add = "jitter") + stat_compare_means() + ylab("B cells") + xlab("")
ggboxplot(temp, "meta", "B_and_plasma_cells", add = "jitter") + stat_compare_means() + ylab("B cells (incl. plasma cells)") + xlab("")

adjustedHRs = data.frame()

for(j in 1:length(kohorty)){
  
  kohorta = kohorty[j]
  temp = all[all$meta == kohorta,]
  #wyniki = list()
  wynikidf = data.frame()
  
  for(i in 1:length(immuoparameters)){
    try({
      selected = immuoparameters[i]
      tempxx = dplyr::select(temp, OS_TIME, OS, Age_at_diagnosis, T, M, Platin, all_of(selected))
      tempxx = tempxx[complete.cases(tempxx),]
      
      f = as.formula(paste0("Surv(OS_TIME, OS) ~ Age_at_diagnosis + 
    T + M + ", selected))
      m = coxph(f , data = tempxx)
      summary(m)
      
      wynikidf[i,"parameter"] = selected
      wynikidf[i,"loglik"] = logLik(m)
      wynikidf[i,"HR"] = summary(m)$conf.int[4,1]
      wynikidf[i,"LowerCI"] = summary(m)$conf.int[4,3]
      wynikidf[i,"UpperCI"] = summary(m)$conf.int[4,4]
      wynikidf[i,"cindex"] = summary(m)$concordance["C"]
      wynikidf[i,"pval"] = summary(m)$coefficients[4,5]
      wynikidf[i,"cohort"] = kohorta
      wynikidf[i,"n"] = sum(temp$Cohort == kohorta)
    })
    
    
    #a = tbl_regression(m, exponentiate = T)
    #wyniki = c(wyniki, list(a))
  }
  
  adjustedHRs = rbind.fill(adjustedHRs, wynikidf)
}
data.table::fwrite(adjustedHRs, "adjustedHRs.csv")

library(meta)
library(metafor)

wyniki_metaanaliz = data.frame()
for(i in 1:length(immuoparameters)) {
  a = try({
    parametr = immuoparameters[i]
    temp = dplyr::filter(adjustedHRs, parameter == parametr)
    temp$logHR = log(temp$HR)
    temp$logUB = log(temp$UpperCI)
    temp$logLB = log(temp$LowerCI)
    temp$selogHR =(temp$logUB - temp$logLB)/(2 * 1.96) 
    
    temp = temp[complete.cases(temp),]
    # temp
    
    metaanalysis <- metagen(logHR, selogHR, studlab = cohort, data = temp, sm = "HR")
    
    
    # groups
    temp$subgroups = unlist(lapply(strsplit(temp$cohort, split = " \\+ "), function(x) x[2]))
    subgroups = unlist(lapply(strsplit(temp$cohort, split = " \\+ "), function(x) x[2]))
    metaregression = metareg(metaanalysis, ~subgroups)
    metaanalysis = update(metaanalysis, byvar = temp$subgroups)
    
    sink(paste0("reports/",parametr,".txt"))
    print(summary(metaanalysis))
    print(summary(metaregression))
    sink()
    
    # summary(metaanalysis)
    # meta::forest(metaanalysis, JAMA.pval = T)
    
    wyniki_metaanaliz[i,"parameter"] = parametr
    wyniki_metaanaliz[i,"HR (fixed)"] =  exp(metaanalysis$TE.fixed)
    wyniki_metaanaliz[i,"Lower CI (fixed)"] =  exp(metaanalysis$lower.fixed)
    wyniki_metaanaliz[i,"Upper CI (fixed)"] =  exp(metaanalysis$upper.fixed)
    wyniki_metaanaliz[i,"p (fixed)"] =  metaanalysis$pval.fixed
    wyniki_metaanaliz[i,"I2"] =  metaanalysis$I2
    wyniki_metaanaliz[i,"HR (random)"] =  exp(metaanalysis$TE.random)
    wyniki_metaanaliz[i,"Lower CI (random)"] =  exp(metaanalysis$lower.random)
    wyniki_metaanaliz[i,"Upper CI (random)"] =  exp(metaanalysis$upper.random)
    wyniki_metaanaliz[i,"p (random)"] =  metaanalysis$pval.random
    
    wyniki_metaanaliz[i,"Common between p"] =  metaanalysis$pval.Q.b.common
    wyniki_metaanaliz[i,"Random between p"] =  metaanalysis$pval.Q.b.random
    
    wyniki_metaanaliz[i,"Metaregression QM"] =  metaregression$QM
    wyniki_metaanaliz[i,"Metaregression p"] =  metaregression$QMp
    
    pdf(paste0("forests/",parametr,".pdf"), width = 20, height = 20)
    meta::forest(metaanalysis)
    #meta::forest(metaregression)
    bubble(metaregression, studlab = T)
    
    dev.off()
    
    # a =a$pval metabias(metaanalysis, method="rank", k.min=6)
    # wyniki_metaanaliz[i,"BeggMazumdar p"] = a$pval
    # 
    # a = metabias(metaanalysis, method="linreg", k.min=6)
    # wyniki_metaanaliz[i,"Egger p"] = 
  })
  if(class(a) == "try-error"){
    file.remove(paste0("forests/",parametr,".pdf"))
  }
}
data.table::fwrite(wyniki_metaanaliz, "metaanalyses.csv")

# CE ----
table(merged_d$B.cells.memory_groups, merged_d$`Carcinoma Ecotype`)
chisq.test(table(merged_d$B.cells.memory_groups, merged_d$`Carcinoma Ecotype`))


temp = dplyr::select(merged_d, starts_with("CE"), B.cells.memory, B_cells)

M = cor(temp)
testRes = cor.mtest(temp, conf.level = 0.95)
library(corrplot)

corrplot(M, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'AOE', diag=FALSE)

OmicSelector::OmicSelector_correlation_plot(temp$B_cells, temp$CE10, "B cells", "CE10", "B cells vs. CE10", yx = F)

# Custom plots ----
library(ggpubr)
library(readr)
merged_d <- read_csv("Files/merged_d.csv")

table(merged_d$Cohort)
library(ggplot2)
library(ggpubr)

library(ggplot2)
library(ggpubr)

# Custom color palette for the Cohorts based on provided image
custom_colors <- c("GREEK" = "#F8766D", "HM" = "#00BA38", "TABER" = "#619CFF")

# Boxplot for TILs
png("Plots/Boxplot_Cohort_TILs.png", width = 6, height = 4, units = "in", res = 300)
ggboxplot(merged_d, x = "Cohort", y = "TILs", add = "jitter", fill = "Cohort") +
  scale_fill_manual(values = custom_colors) +
  geom_pwc(tip.length = 0, label = "p = {p.adj.format}", method = "t_test") +
  theme_minimal()
dev.off()

# Boxplot for B cells
png("Plots/Boxplot_Cohort_Bcells.png", width = 6, height = 4, units = "in", res = 300)
ggboxplot(merged_d, x = "Cohort", y = "B_cells", add = "jitter", fill = "Cohort") +
  scale_fill_manual(values = custom_colors) +
  geom_pwc(tip.length = 0, label = "p = {p.adj.format}", method = "t_test") +
  ylab("B cells") +
  theme_minimal()
dev.off()

# Boxplot for Lymphoid lineage
png("Plots/Boxplot_Cohort_Lymphoid.png", width = 6, height = 4, units = "in", res = 300)
ggboxplot(merged_d, x = "Cohort", y = "Lymphoid", add = "jitter", fill = "Cohort") +
  scale_fill_manual(values = custom_colors) +
  geom_pwc(tip.length = 0, label = "p = {p.adj.format}", method = "t_test") +
  ylab("Lymphoid lineage") +
  theme_minimal()
dev.off()

# Boxplot for B cells memory
png("Plots/Boxplot_Cohort_BcellsMem.png", width = 6, height = 4, units = "in", res = 300)
ggboxplot(merged_d, x = "Cohort", y = "B.cells.memory", add = "jitter", fill = "Cohort") +
  scale_fill_manual(values = custom_colors) +
  geom_pwc(tip.length = 0, label = "p = {p.adj.format}", method = "t_test") +
  ylab("B cells memory") +
  theme_minimal()
dev.off()


library(OmicSelector)
png("Plots/Correlation_CE10_BcellsMem.png", width = 6, height = 5, units = "in", res = 300)
OmicSelector_correlation_plot(merged_d$CE10, merged_d$B.cells.memory,"CE10 score", "B cells memory",  "", yx=F)
dev.off()

png("Plots/Correlation_CE6_BcellsMem.png", width = 6, height = 5, units = "in", res = 300)
OmicSelector_correlation_plot(merged_d$CE6, merged_d$B.cells.memory,"CE6 score", "B cells memory",  "", yx=F)
dev.off()

png("Plots/Correlation_CE7_BcellsMem.png", width = 6, height = 5, units = "in", res = 300)
OmicSelector_correlation_plot(merged_d$CE7, merged_d$B.cells.memory,"CE7 score", "B cells memory",  "", yx=F)
dev.off()

png("Plots/Correlation_LumNS_BcellsMem.png", width = 6, height = 5, units = "in", res = 300)
OmicSelector_correlation_plot(merged_d$LumNS, merged_d$B.cells.memory,"LumNS", "B cells memory",  "", yx=F, metoda = "spearman")
dev.off()

png("Plots/Correlation_Neutrophils_BcellsMem.png", width = 6, height = 5, units = "in", res = 300)
OmicSelector_correlation_plot(merged_d$Neutrophils, merged_d$B.cells.memory,"Neutrophils", "B cells memory",  "", yx=F)
dev.off()

merged_d = data.table::fread("Files/merged_d.csv")


# Improv 10.10.2024 ----
setwd("H:/My Drive/PostCDDP/ANALYSIS/v4_032024")
all = data.table::fread("SelectedSignatures/all.csv.gz")

# Load necessary libraries
library(ggplot2)
library(ggpubr)

all$B.cells.memory

# Scatter plot with correlation and linear regression line
plot1 <- ggplot(all, aes(x = bcell_signature1, y = B.cells.memory)) +
  geom_point() +
  geom_smooth(method = "lm", color = "darkred", fill = "gray") +
  stat_cor(method = "spearman", label.x = 0, label.y = 1, p.accuracy = 0.001, r.accuracy = 0.01) + # Adjust label positions as needed
  labs(title = "",
       x = "IFNG Expression / Total CD8 Fraction",
       y = "B-cell memory (CIBERSORTx)") + ylim(0,1) +
  theme_bw()
plot1

plot2 <- ggplot(data, aes(x = IFNG_CD8_Fraction, y = ICB_Exposed_TAM_Sig_Score)) +
  geom_point() +
  geom_smooth(method = "lm", color = "purple", fill = "gray") +
  stat_cor(method = "pearson", label.x = 5, label.y = 2.8, digits = 4) + # Adjust label positions as needed
  labs(x = "IFNG Expression / Total CD8 Fraction",
       y = "ICB-Exposed TAM Sig. Score") +
  theme_minimal()

# Arrange the plots side by side
combined_plot <- ggarrange(plot1, plot2, ncol = 2, nrow = 1)
print(combined_plot)
