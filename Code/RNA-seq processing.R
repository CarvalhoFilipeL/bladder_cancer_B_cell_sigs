####################################################
# 1. Setup & Load Libraries
####################################################
setwd("/home/kgs24/PostCDDP/Analysis_NF/DiffExp")

# Core libraries for data import and processing
library(readxl)
library(tidyverse)
library(data.table)
library(tximport)      # Import transcript-level quantification
library(readr)
library(DESeq2)        # Differential expression and filtering
library(edgeR)         # TMM normalization
library(consensusMIBC) # Consensus mRNA clustering
library(BLCAsubtyping) # mRNA subtyping
library(immunedeconv)  # Immune deconvolution (here: CIBERSORTx)
library(Biobase)       # For building ExpressionSet objects

####################################################
# 2. Read Metadata & Define Sample Information
####################################################
# Read metadata containing sample names and other info.
metadata <- fread("../metasheet.csv")

# Create a sample table with paths to quantification files
cData <- metadata %>%
  mutate(
    names = sample,
    files = paste0("../Quantification/star_salmon/", sample, "/quant.sf"),
    Batch = as.factor(ifelse(startsWith(sample, "hm"), "HM", "GREEK"))
  )
table(cData$Batch)

####################################################
# 3. Import Quantification Data with tximport
####################################################
# Read transcript-to-gene mapping file.
tx_ref <- read_table("../Quantification/star_salmon/salmon_tx2gene.tsv", 
                     col_names = FALSE)
colnames(tx_ref) <- c("TRANSCRIPT_ID", "GENE_ID", "GENE_SYMBOL")
head(tx_ref)

# Import quantifications (Salmon outputs) and summarize to gene level.
txi <- tximport(cData$files,
                type = "salmon",
                tx2gene = tx_ref[, 1:2],
                countsFromAbundance = "lengthScaledTPM")

####################################################
# 4. Create DESeq2 Dataset & Filter Low-Count Genes
####################################################
# Build DESeq2 dataset using the imported counts.
dds <- DESeqDataSetFromTximport(txi = txi, colData = cData, design = ~ Batch)

# Filter out genes with a total read count less than 10.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

####################################################
# 5. Normalize with edgeR (TMM Normalization)
####################################################
# Create a DGEList object from raw counts.
y <- DGEList(counts = counts(dds), group = rep("Unknown", ncol(dds)))
y <- calcNormFactors(y)

# Compute log2-transformed counts-per-million (CPM) using TMM normalization.
TMM <- cpm(y, normalized.lib.sizes = TRUE, log = TRUE)
dim(TMM)

# Save TMM-normalized data and raw counts.
colnames(TMM) <- metadata$sample
write.csv(data.frame("Gene" = rownames(TMM), TMM),
          "Data/TMM.csv", quote = FALSE, row.names = FALSE)
# For methods (e.g., CIBERSORTx) requiring linear-scale input, square TMM.
write.table(data.frame("Gene" = rownames(TMM), TMM^2),
            "Data/TMM_for_cibersort.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
write.csv(counts(dds), "Data/counts.csv", quote = FALSE)

####################################################
# 6. Consensus mRNA Clustering & Subtyping
####################################################
# Generate consensus mRNA clusters using the normalized TMM data.
if (!dir.exists("mRNA_Clusters")) { dir.create("mRNA_Clusters") }
clusters <- getConsensusClass(TMM, minCor = 0.2, gene_id = "hgnc_symbol")
fwrite(clusters, "mRNA_Clusters/clusters.csv")
table(clusters$consensusClass)

# Run additional mRNA subtyping using BLCAsubtyping.
cl <- classify(expMat = TMM,
               classification.systems = c("Baylor", "UNC", "MDA", "CIT", "Lund", "TCGA"))
fwrite(cl, "mRNA_Clusters/subtypes.csv")
fwrite(cl, "Summarized/subtypes.csv")

####################################################
# 7. Immune Infiltration Using CIBERSORTx
####################################################
# For immune deconvolution, we use only CIBERSORTx.
# CIBERSORTx requires linear-scale data. Here, we use TMM^2.
x <- TMM^2
# Ensure sample names are valid.
colnames(x) <- make.names(colnames(x), unique = TRUE)

# (If required, set the location of the CIBERSORTx binary and LM22 matrix.
#  For example:
# set_cibersort_binary("/home/kgs24/PostCDDP/Soft/CIBERSORT/CIBERSORT.R")
# set_cibersort_mat("/home/kgs24/PostCDDP/Soft/CIBERSORT/LM22.txt")
# )

# Run CIBERSORTx (absolute mode is used here; adjust the 'method' if needed)
method <- "cibersort_abs"
cibersortx_res <- deconvolute(x, method)
# Transpose and clean up the output so that the first row becomes column names.
cibersortx_res <- t(cibersortx_res)
colnames(cibersortx_res) <- cibersortx_res[1,]
cibersortx_res <- cibersortx_res[-1,]
cibersortx_res <- as_tibble(cibersortx_res)
colnames(cibersortx_res) <- paste0("CIBERSORTx_", colnames(cibersortx_res))
# Save CIBERSORTx results.
fwrite(cibersortx_res, "Immune_infilitration/CIBERSORTx.csv", col.names = TRUE)

####################################################
# 8. Merge & Summarize Metadata
####################################################
# Merge various summary files (e.g., mRNA clusters, subtypes, and immune infiltration)
# with the main metadata for downstream analyses.
metadane <- read_excel("/home/kgs24/PostCDDP/v1_Dec2022/PostCDDP_all_metadata.xlsx", 
                       sheet = "HM+GREEK")

# Merge in consensus mRNA clusters and subtypes.
temp <- fread("mRNA_Clusters/clusters.csv")
metadane <- cbind(metadane, temp)
temp <- fread("mRNA_Clusters/subtypes.csv")
metadane <- cbind(metadane, temp)

# Merge in immune deconvolution results (CIBERSORTx).
temp <- fread("Immune_infilitration/CIBERSORTx.csv")
metadane <- cbind(metadane, temp)

# Remove redundant columns (if any) and duplicate entries.
metadane <- metadane[!duplicated(as.list(metadane))]
metadane <- select(metadane, -starts_with("zz"))
colnames(metadane) <- paste0("meta_", colnames(metadane))
fwrite(metadane, "metadata.csv")

# Save the workspace image for checkpointing.
save.image("checkpoint.RData")
