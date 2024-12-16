
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  install.packages("AnnotationDbi")
}
if (!requireNamespace("hgu133plus2.db", quietly = TRUE)) {
  BiocManager::install("hgu133plus2.db")
}
install.packages("devtools")
devtools::install_github("kevinblighe/CorLevelPlot",force = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  BiocManager::install("WGCNA",force=TRUE)
}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("impute")
library(impute)
library('WGCNA')
library(biomaRt)
install.packages("doParallel")
install.packages("foreach")
library(doParallel)
library(foreach)
install.packages("BiocManager")
 BiocManager::install("beadarray")
library(CorLevelPlot)
 
# Load required libraries
library(AnnotationDbi)
library(hgu133plus2.db)
library(CorLevelPlot)
BiocManager::install("KEGGREST")
library(KEGGREST)
library(KEGG.db)
library(affy)
library(oligo)
library(Biobase)
library(GEOquery)
library(dplyr)
library(clusterProfiler)
BiocManager::install("clusterProfiler",force = TRUE)
BiocManager::install("arrayQualityMetrics")
BiocManager::install("VennDiagram")
library(arrayQualityMetrics)
install.packages("splitstackshape")
library("splitstackshape")
library("tidyr")
library(limma)
library('org.Hs.eg.db')
library(dplyr)
library(VennDiagram)

if (!requireNamespace("hugene10sttranscriptcluster.db", quietly = TRUE)) {
  BiocManager::install("hugene10sttranscriptcluster.db")
}
library(hugene10sttranscriptcluster.db)

#untar the file
untar("GSE32537_RAW.tar")
celFiles <- list.celfiles()
affyRAW <- read.celfiles(celFiles)
hist(affyRAW, target = "probeset")
eset <- oligo::rma(affyRAW)
hist(eset)
write.exprs(eset,file="data1.txt")

#annotation
mydata <- read.delim("GSE32537_family.soft", check.names = FALSE)
abc <- data.frame(mydata)
abc$Gene_Symbol <- sapply(strsplit(as.character(abc$gene_assignment), " // "), function(x) trimws(x[2]))

#removing NA values from GeneSymbol column
abc <- abc[!is.na(abc$Gene_Symbol), ]
extracted_col=abc$Gene_Symbol
abc$Gene_Symbol <-abc$Gene_Symbol
write.table(abc, file="imp1.txt",sep = "\t", row.names=FALSE,quote = FALSE)
data1 <- read.delim("imp1.txt",check.names = FALSE)
data2 <- read.delim("data1.txt", check.names=FALSE)

#to check the columns in both data
colnames(data1)
colnames(data2)

# Convert both ID columns to character
data1$ID <- as.character(data1$ID)
data2$ID <- as.character(data2$ID)


# Perform the join
combined <- left_join(data1, data2, by = "ID")

write.csv(combined,"annotated.csv")
pData(eset)
targets <- readTargets("target_modified.csv", sep = ",")

# Set "control" as the reference level and create the design matrix
unique(targets$condition)

# Define conditions with valid column names
condition <- factor(targets$condition, levels = c("control", "COP", "RB_ILD", "DIP", "NSIP", "UF", "IPF_UIP"))
design <- model.matrix(~ 0 + condition)

# Modify column names to valid R variable names
colnames(design) <- c("control", "COP", "RB_ILD", "DIP", "NSIP", "UF", "IPF_UIP")
print(design)

rownames(design)
colnames(eset)

dim(design) # Should be (samples, factors)
dim(eset)   # Should be (genes, samples)


# Fit linear model
fit <- lmFit(eset, design)

# Differential expression analysis for COP vs control
contrast.matrix_COP <- makeContrasts(COP - control, levels = design)
fit_COP <- contrasts.fit(fit, contrast.matrix_COP)
fit_COP <- eBayes(fit_COP)
res_COP_vs_control <- topTable(fit_COP, number=Inf)
write.table(res_COP_vs_control, "COP_vs_control.txt", sep="\t")

# Differential expression analysis for RB_ILD vs control
contrast.matrix_RB_ILD <- makeContrasts(RB_ILD - control, levels = design)
fit_RB_ILD <- contrasts.fit(fit, contrast.matrix_RB_ILD)
fit_RB_ILD <- eBayes(fit_RB_ILD)
res_RB_ILD_vs_control <- topTable(fit_RB_ILD, number=Inf)
write.table(res_RB_ILD_vs_control, "RB_ILD_vs_control.txt", sep="\t")

# Differential expression analysis for DIP vs control
contrast.matrix_DIP <- makeContrasts(DIP - control, levels = design)
fit_DIP <- contrasts.fit(fit, contrast.matrix_DIP)
fit_DIP <- eBayes(fit_DIP)
res_DIP_vs_control <- topTable(fit_DIP, number=Inf)
write.table(res_DIP_vs_control, "DIP_vs_control.txt", sep="\t")

# Differential expression analysis for NSIP vs control
contrast.matrix_NSIP <- makeContrasts(NSIP - control, levels = design)
fit_NSIP <- contrasts.fit(fit, contrast.matrix_NSIP)
fit_NSIP <- eBayes(fit_NSIP)
res_NSIP_vs_control <- topTable(fit_NSIP, number=Inf)
write.table(res_NSIP_vs_control, "NSIP_vs_control.txt", sep="\t")

# Differential expression analysis for UF vs control
contrast.matrix_UF <- makeContrasts(UF - control, levels = design)
fit_UF <- contrasts.fit(fit, contrast.matrix_UF)
fit_UF <- eBayes(fit_UF)
res_UF_vs_control <- topTable(fit_UF, number=Inf)
write.table(res_UF_vs_control, "UF_vs_control.txt", sep="\t")

# Differential expression analysis for IPF_UIP vs control
contrast.matrix_IPF_UIP <- makeContrasts(IPF_UIP - control, levels = design)
fit_IPF_UIP <- contrasts.fit(fit, contrast.matrix_IPF_UIP)
fit_IPF_UIP <- eBayes(fit_IPF_UIP)
res_IPF_UIP_vs_control <- topTable(fit_IPF_UIP, number=Inf)
write.table(res_IPF_UIP_vs_control, "IPF_UIP_vs_control.txt", sep="\t")

# Output results are saved to 'COP_vs_control.txt', 'DIP_vs_control.txt','NSIP_vs_control.txt', 'RB_ILD_vs_control.txt', 'UF_vs_control.txt', 'IPF_UIP_vs_control.txt'

# Load your DE results from the three comparisons
res1 <- read.delim("COP_vs_control.txt", check.names = FALSE)
res2 <- read.delim("DIP_vs_control.txt", check.names = FALSE)
res3 <- read.delim("NSIP_vs_control.txt", check.names = FALSE)
res4 <- read.delim("RB_ILD_vs_control.txt", check.names = FALSE)
res5 <- read.delim("UF_vs_control.txt", check.names = FALSE)
res6 <- read.delim("IPF_UIP_vs_control.txt", check.names = FALSE)


# Define significance threshold
significance_threshold <- 0.05
logFC_threshold <- 1

# Filter DE genes for each dataset
de_genes1 <- rownames(res1[res1$P.Value < significance_threshold & abs(res1$logFC) > logFC_threshold, ])
de_genes2 <- rownames(res2[res2$P.Value < significance_threshold & abs(res2$logFC) > logFC_threshold, ])
de_genes3 <- rownames(res3[res3$P.Value < significance_threshold & abs(res3$logFC) > logFC_threshold, ])
de_genes4 <- rownames(res4[res4$P.Value < significance_threshold & abs(res4$logFC) > logFC_threshold, ])
de_genes5 <- rownames(res5[res5$P.Value < significance_threshold & abs(res5$logFC) > logFC_threshold, ])
de_genes6 <- rownames(res6[res6$P.Value < significance_threshold & abs(res6$logFC) > logFC_threshold, ])

# Convert to character vectors
de_genes1 <- as.character(de_genes1)
de_genes2 <- as.character(de_genes2)
de_genes3 <- as.character(de_genes3)
de_genes4 <- as.character(de_genes4)
de_genes5 <- as.character(de_genes5)
de_genes6 <- as.character(de_genes6)

# Find common DE genes
common_de_genes <- Reduce(intersect, list(de_genes1, de_genes2, de_genes3,de_genes4,de_genes5,de_genes6))

# Print common DE genes
print("Common Differentially Expressed Genes:")
print(common_de_genes)

##########################################################
#pathway enrichment

# Load the differential expression results (from earlier analysis)
res_COP_vs_control<- read.delim("COP_vs_control.txt", check.names = FALSE)
res_RB_ILD_vs_control <- read.delim("RB_ILD_vs_control.txt", check.names = FALSE)
res_DIP_vs_control <- read.delim("DIP_vs_control.txt", check.names = FALSE)
res_NSIP_vs_control <- read.delim("NSIP_vs_control.txt", check.names = FALSE)
res_UF_vs_control <- read.delim("UF_vs_control.txt", check.names = FALSE)
res_IPF_UIP_vs_control <- read.delim("IPF_UIP_vs_control.txt", check.names = FALSE)

# Annotate DEGs according to logFC and adjusted P-value
# For COP vs control
res_COP_vs_control <- res_COP_vs_control %>%
  mutate(diffexpressed = case_when(
    logFC > 0 & adj.P.Val < 0.05 ~ 'UP',
    logFC < 0 & adj.P.Val < 0.05 ~ 'DOWN',
    adj.P.Val > 0.05 ~ 'NO'
  ))

# For RB_ILD_vs_control
res_RB_ILD_vs_control <- res_RB_ILD_vs_control %>%
  mutate(diffexpressed = case_when(
    logFC > 0 & adj.P.Val < 0.05 ~ 'UP',
    logFC < 0 & adj.P.Val < 0.05 ~ 'DOWN',
    adj.P.Val > 0.05 ~ 'NO'
  ))

# For DIP_vs_control
res_DIP_vs_control <- res_DIP_vs_control %>%
  mutate(diffexpressed = case_when(
    logFC > 0 & adj.P.Val < 0.05 ~ 'UP',
    logFC < 0 & adj.P.Val < 0.05 ~ 'DOWN',
    adj.P.Val > 0.05 ~ 'NO'
  ))
# For NSIP_vs_control
res_NSIP_vs_control <- res_NSIP_vs_control %>%
  mutate(diffexpressed = case_when(
    logFC > 0 & adj.P.Val < 0.05 ~ 'UP',
    logFC < 0 & adj.P.Val < 0.05 ~ 'DOWN',
    adj.P.Val > 0.05 ~ 'NO'
  ))
# For UF_vs_control
res_UF_vs_control <- res_UF_vs_control %>%
  mutate(diffexpressed = case_when(
    logFC > 0 & adj.P.Val < 0.05 ~ 'UP',
    logFC < 0 & adj.P.Val < 0.05 ~ 'DOWN',
    adj.P.Val > 0.05 ~ 'NO'
  ))
# For IPF_UIP_vs_control
res_IPF_UIP_vs_control <- res_IPF_UIP_vs_control %>%
  mutate(diffexpressed = case_when(
    logFC > 0 & adj.P.Val < 0.05 ~ 'UP',
    logFC < 0 & adj.P.Val < 0.05 ~ 'DOWN',
    adj.P.Val > 0.05 ~ 'NO'
  ))


# Now save the annotated files if you want
write.csv(res_COP_vs_control, "COP_vs_control_annotated.csv")
write.csv(res_RB_ILD_vs_control, "RB_ILD_vs_control_annotated.csv")
write.csv(res_DIP_vs_control, "DIP_vs_control_annotated.csv")
write.csv(res_NSIP_vs_control, "NSIP_vs_control_annotated.csv")
write.csv(res_UF_vs_control, "UF_vs_control_annotated.csv")
write.csv(res_IPF_UIP_vs_control, "IPF_UIP_vs_control_annotated.csv")

#########################
# Load the annotated CSV files
COP_vs_control_annotated <- read.csv("COP_vs_control_annotated.csv", row.names = 1)
RB_ILD_vs_control_annotated <- read.csv("RB_ILD_vs_control_annotated.csv", row.names = 1)
DIP_vs_control_annotated <- read.csv("DIP_vs_control_annotated.csv", row.names = 1)
NSIP_vs_control_annotated<- read.csv("NSIP_vs_control_annotated.csv", row.names = 1)
UF_vs_control_annotated <- read.csv("UF_vs_control_annotated.csv", row.names = 1)
IPF_UIP_vs_control_annotated <- read.csv("IPF_UIP_vs_control_annotated.csv", row.names = 1)

####################################################

# For COP vs. control
probe_ids_COP <- rownames(res_COP_vs_control)
gene_symbols_COP <- mapIds(hugene10sttranscriptcluster.db, keys = probe_ids_COP, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
res_COP_vs_control$gene_symbol <- gene_symbols_COP

# For RB_ILD_vs_control
probe_ids_RB_ILD <- rownames(res_RB_ILD_vs_control)
gene_symbols_RB_ILD <- mapIds(hugene10sttranscriptcluster.db, keys = probe_ids_RB_ILD, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
res_RB_ILD_vs_control$gene_symbol <- gene_symbols_RB_ILD

# For DIP_vs_control
probe_ids_DIP <- rownames(res_DIP_vs_control)
gene_symbols_DIP <- mapIds(hugene10sttranscriptcluster.db, keys = probe_ids_DIP, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
res_DIP_vs_control$gene_symbol <- gene_symbols_DIP

# For NSIP_vs_control
probe_ids_NSIP <- rownames(res_NSIP_vs_control)
gene_symbols_NSIP <- mapIds(hugene10sttranscriptcluster.db, keys = probe_ids_NSIP, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
res_NSIP_vs_control$gene_symbol <- gene_symbols_NSIP

# For UF_vs_control
probe_ids_UF <- rownames(res_UF_vs_control)
gene_symbols_UF <- mapIds(hugene10sttranscriptcluster.db, keys = probe_ids_UF, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
res_UF_vs_control$gene_symbol <- gene_symbols_UF

# For IPF_UIP_vs_control
probe_ids_IPF_UIP <- rownames(res_IPF_UIP_vs_control)
gene_symbols_IPF_UIP <- mapIds(hugene10sttranscriptcluster.db, keys = probe_ids_IPF_UIP, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
res_IPF_UIP_vs_control$gene_symbol <- gene_symbols_IPF_UIP

# Extract probe IDs
probe_ids <- rownames(res_COP_vs_control)

# Get gene symbols
gene_symbols_COP <- mapIds(hugene10sttranscriptcluster.db, keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")

# Add gene symbols to your data frame
res_COP_vs_control$gene_symbol <- gene_symbols_COP

# View the updated data frame with gene symbols
head(res_COP_vs_control)
head(res_RB_ILD_vs_control)
head(res_DIP_vs_control)
head(res_NSIP_vs_control)
head(res_UF_vs_control)
head(res_IPF_UIP_vs_control)

##########################################################

common_de_genes
### common_de_genes had only probe ids Gotta convert it into entrez first


# Convert probe IDs to gene symbols
gene_symbols <- mapIds(hugene10sttranscriptcluster.db, keys = common_de_genes, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")

entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Check the Entrez IDs
print(entrez_ids)

# KEGG Pathway Enrichment Analysis (using Entrez IDs)
kegg_enrich <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = "hsa")
kegg_results_df <- as.data.frame(kegg_enrich)
print(kegg_results_df)



# GO Enrichment Analysis (Biological Process)
go_enrich <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", 
                      ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)


# Convert GO enrichment results to a data frame
go_enrich_df <- as.data.frame(go_enrich)

# View the first few rows of the results
head(go_enrich_df)


common_de_genes
gene_symbols <- na.omit(gene_symbols)
##################################################
#visualization
# Assuming 'common_de_genes' is your gene list
# Use org.Hs.eg.db for human genes; change for other organisms
go_enrich <- enrichGO(gene = gene_symbols,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",  # Biological Process (can also use "MF" or "CC")
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

# View top GO terms (optional) illitanka aaju ishappp kunniiiii
head(go_enrich)

gene_symbols

barplot(go_enrich, showCategory = 10, title = "Top 10 GO Terms Enriched")
dotplot(go_enrich, showCategory = 10, title = "Top 10 GO Terms Enriched")

cnetplot(go_enrich, showCategory = 5, foldChange = NULL, layout = "kk", title = "GO Enrichment Map")


#kegg enrichment 
# Perform KEGG enrichment analysis using Entrez IDs

entrez_ids <- na.omit(entrez_ids)


# Perform KEGG enrichment analysis using the enricher() function
entrez_ids

kk <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)


# Barplot for the top 10 enriched KEGG pathways
barplot(kk, showCategory = 10, title = "Top 10 KEGG Pathways Enriched")

# Dotplot for enriched pathways
dotplot(kk, showCategory = 10, title = "Top 10 KEGG Pathways Enriched")

# this section is GSEA
# Ranked list based on logFC for COP vs control
geneList_COP <- res_COP_vs_control$logFC
names(geneList_COP) <- res_COP_vs_control$gene_symbol
geneList_COP <- sort(geneList_COP, decreasing = TRUE)

# Ranked list for RB_ILD_vs_control
geneList_RB_ILD <- res_RB_ILD_vs_control$logFC
names(geneList_RB_ILD) <- res_RB_ILD_vs_control$gene_symbol
geneList_RB_ILD <- sort(geneList_RB_ILD, decreasing = TRUE)

# Ranked list for DIP_vs_control
geneList_DIP <- res_DIP_vs_control$logFC
names(geneList_DIP) <- res_DIP_vs_control$gene_symbol
geneList_DIP <- sort(geneList_DIP, decreasing = TRUE)

# Ranked list for NSIP_vs_control
geneList_NSIP <- res_NSIP_vs_control$logFC
names(geneList_NSIP) <- res_NSIP_vs_control$gene_symbol
geneList_NSIP <- sort(geneList_NSIP, decreasing = TRUE)

# Ranked list for UF_vs_control
geneList_UF <- res_UF_vs_control$logFC
names(geneList_UF) <- res_UF_vs_control$gene_symbol
geneList_UF <- sort(geneList_UF, decreasing = TRUE)

# Ranked list for IPF_UIP_vs_control
geneList_IPF_UIP <- res_IPF_UIP_vs_control$logFC
names(geneList_IPF_UIP) <- res_IPF_UIP_vs_control$gene_symbol
geneList_IPF_UIP <- sort(geneList_IPF_UIP, decreasing = TRUE)


# GSEA for COP vs control (GO terms)
gsea_GO_COP <- gseGO(geneList = geneList_COP, 
                         OrgDb = org.Hs.eg.db, 
                         keyType = "SYMBOL", 
                         ont = "BP", 
                         pvalueCutoff = 0.05, 
                         verbose = FALSE)
head(gsea_GO_COP)
# GSEA for RB_ILD_vs_control (GO terms)
gsea_GO_RB_ILD <- gseGO(geneList = geneList_RB_ILD, 
                     OrgDb = org.Hs.eg.db, 
                     keyType = "SYMBOL", 
                     ont = "BP", 
                     pvalueCutoff = 0.05, 
                     verbose = FALSE)

# GSEA for DIP_vs_control (GO terms)
gsea_GO_DIP <- gseGO(geneList = geneList_DIP, 
                     OrgDb = org.Hs.eg.db, 
                     keyType = "SYMBOL", 
                     ont = "BP", 
                     pvalueCutoff = 0.05, 
                     verbose = FALSE)
# GSEA for NSIP_vs_control (GO terms)
gsea_GO_NSIP <- gseGO(geneList = geneList_NSIP, 
                     OrgDb = org.Hs.eg.db, 
                     keyType = "SYMBOL", 
                     ont = "BP", 
                     pvalueCutoff = 0.05, 
                     verbose = FALSE)

# GSEA for UF_vs_control (GO terms)
gsea_GO_UF <- gseGO(geneList = geneList_UF, 
                     OrgDb = org.Hs.eg.db, 
                     keyType = "SYMBOL", 
                     ont = "BP", 
                     pvalueCutoff = 0.05, 
                     verbose = FALSE)

geneList_IPF_UIP <- c("IQCN" = 3.625597, "MRE11" = 3.121958, "UCHL1" = 3.071618, "TDRD9" = 2.832898, "MLH3" = 2.803998, "SKIL" = 2.688647)
mapped_genes <- mapIds(org.Hs.eg.db, keys = names(geneList_IPF_UIP), column = "SYMBOL", keytype = "SYMBOL")
names(geneList_IPF_UIP)

# GSEA for IPF_UIP_vs_control (GO terms)
gsea_GO_IPF_UIP <- gseGO(geneList = geneList_IPF_UIP, 
                     OrgDb = org.Hs.eg.db, 
                     keyType = "SYMBOL", 
                     ont = "BP", 
                     pvalueCutoff = 0.05, 
                     verbose = FALSE)

head(geneList_IPF_UIP)



dotplot(gsea_GO_COP, showCategory = 10, title = "GSEA GO COP vs control")
 

#RB_ILD vs control

dotplot(gsea_GO_RB_ILD, showCategory = 10, title = "GSEA GO RB_ILD vs control")

# DIP_vs_control

dotplot(gsea_GO_DIP, showCategory = 10, title = "GSEA GO DIP vs control")

#NSIP_vs_control
dotplot(gsea_GO_NSIP, showCategory = 10, title = "GSEA GO NSIP vs control")
#UF_vs_control
dotplot(gsea_GO_UF, showCategory = 10, title = "GSEA GO UF vs control")

#IPF_UIP_vs_control
dotplot(gsea_GO_IPF_UIP, showCategory = 10, title = "GSEA GO IPF_UIP vs control")


############################################################################
#WGCNA Analysis

############################################################################
options(stringsAsFactors = FALSE)
expression_data <- as.data.frame(exprs(eset))  # Extract expression data from `eset`
expression_data <- t(expression_data)  # Transpose to get genes as columns, samples as rows


gsg = goodSamplesGenes(expression_data)
summary(gsg)
gsg$allOK
table(gsg$goodGenes)


# Calculate variance for each gene
gene_variances <- apply(expression_data, 2, var)

# Set a threshold to keep genes with higher variability (e.g., top 50%)
threshold <- quantile(gene_variances, 0.5)
filtered_expression_data <- expression_data[, gene_variances > threshold]

dim(filtered_expression_data)

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  BiocManager::install("WGCNA")
}
library(WGCNA)


options(stringsAsFactors = FALSE)
expression_data <- as.data.frame(exprs(eset))  # Extract expression data from `eset`
expression_data <- t(expression_data)  # Transpose to get genes as columns, samples as rows


gsg = goodSamplesGenes(expression_data)
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
powers = c(1:20)
sft = pickSoftThreshold(
  filtered_expression_data,
  powerVector = powers,
  verbose = 5,
  networkType = "signed"
)


# Plot Scale-Free Topology Fit
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.8, col = "blue",ity=2)

# Plot Mean Connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "b", pch = 19, col = "blue", main = "Mean Connectivity vs Power")

# Calculate variance for each gene
gene_variances <- apply(expression_data, 2, var)

# Set a threshold to keep genes with higher variability (e.g., top 50%)
threshold <- quantile(gene_variances, 0.5)
filtered_expression_data <- expression_data[, gene_variances > threshold]

dim(filtered_expression_data)

# Convert expression_data to norm.counts
norm.counts <- filtered_expression_data
norm.counts[] = sapply(norm.counts, as.numeric)  # Convert to numeric

soft_power <- 8
temp_cor <- cor
cor <- WGCNA::cor  # Temporarily reassign cor to WGCNA’s version

# Execute blockwiseModules
n_threads <- 16  # Use the full 16 threads

module_colors <- blockwiseModules(
  norm.counts, maxBlockSize = 14000, TOMType = "signed",
  power = soft_power, mergeCutHeight = 0.25,
  numericLabels = FALSE, randomSeed = 1234, verbose = 3,
  nThreads = n_threads  # Enable multithreading
)

cor <- temp_cor

# Get the block assignment for each gene
blocks <- module_colors$blocks
print("Number of blocks:")
print(table(blocks))

# Plot dendrogram for each block separately
for(block in unique(blocks)) {
  # Get genes in this block
  block_genes <- which(blocks == block)
  
  # Get corresponding colors for this block
  block_colors <- data.frame(
    unmerged = module_colors$unmergedColors[block_genes],
    merged = module_colors$colors[block_genes]
  )
  
  # Create plot title
  plot_title <- paste("Block", block, "Dendrogram and Module Colors")
  
  # Plot dendrogram for this block
  pdf(paste0("Block_", block, "_dendrogram.pdf"), width = 12, height = 8)
  plotDendroAndColors(
    module_colors$dendrograms[[block]], 
    block_colors,
    c("unmerged", "merged"),
    main = plot_title,
    dendroLabels = FALSE,
    addGuide = TRUE,
    hang = 0.03,
    guideHang = 0.05
  )
  dev.off()
}

# Get the block assignment for each gene
blocks <- module_colors$blocks
print("Number of blocks:")
print(table(blocks))

# Plot dendrogram for each block separately
par(mfrow = c(length(unique(blocks)), 1))  # Set up plotting layout
for(block in unique(blocks)) {
  # Get genes in this block
  block_genes <- which(blocks == block)
  
  # Get corresponding colors for this block
  block_colors <- data.frame(
    unmerged = module_colors$unmergedColors[block_genes],
    merged = module_colors$colors[block_genes]
  )
  
  # Create plot title
  plot_title <- paste("Block", block, "Dendrogram and Module Colors")
  
  # Plot dendrogram for this block
  plotDendroAndColors(
    module_colors$dendrograms[[block]], 
    block_colors,
    c("unmerged", "merged"),
    main = plot_title,
    dendroLabels = FALSE,
    addGuide = TRUE,
    hang = 0.03,
    guideHang = 0.05
  )
}
par(mfrow = c(1, 1))  # Reset plotting layout

# Module Eigengenes analysis
module_eigengene <- module_colors$MEs
head(module_eigengene)

# Print gene counts for each module
print("Module sizes:")
print(table(module_colors$colors))

# Create diagnostic checks output
cat("Number of genes in dendrogram:", length(module_colors$dendrograms[[1]]$order), "\n")
cat("Number of genes in unmergedColors:", length(module_colors$unmergedColors), "\n")
cat("Number of genes in colors:", length(module_colors$colors), "\n")

# Plot final dendrogram with color matching
colors_to_plot <- data.frame(
  unmerged = module_colors$unmergedColors[module_colors$dendrograms[[1]]$order],
  merged = module_colors$colors[module_colors$dendrograms[[1]]$order]
)

plotDendroAndColors(
  module_colors$dendrograms[[1]], 
  colors_to_plot,
  c("unmerged", "merged"),
  dendroLabels = FALSE,
  addGuide = TRUE,
  hang = 0.03,
  guideHang = 0.05
)

blocks <- module_colors$blocks
print("Number of blocks:")
print(table(blocks))

# Plot dendrogram for each block in R graphics window
for (block in unique(blocks)) {
  block_genes <- which(blocks == block)
  block_colors <- data.frame(
    unmerged = module_colors$unmergedColors[block_genes],
    merged = module_colors$colors[block_genes]
  )
  
  plot_title <- paste("Block", block, "Dendrogram and Module Colors")
  plotDendroAndColors(
    module_colors$dendrograms[[block]], 
    block_colors,
    c("unmerged", "merged"),
    main = plot_title,
    dendroLabels = FALSE,
    addGuide = TRUE,
    hang = 0.03,
    guideHang = 0.05
  )
}

#### Module Eigens 
module_eigengene <- module_colors$MEs


head(module_eigengene)

targets$condition <- factor(targets$condition,levels = c("control", "COP", "RB_ILD", "DIP", "NSIP", "UF", "IPF_UIP"))
levels(targets$condition)
disease_type = binarizeCategoricalColumns(targets$condition,
                           includePairwise=FALSE,
                           includeLevelVsAll = TRUE,
                           minCount = 1)
traits = targets
traits = cbind(traits,disease_type)
colnames(traits)
traits <- traits[, !(names(traits) %in% "condition")]

targets
# defining no. of genes and samples
nSamples = nrow(norm.counts)


nGenes = ncol(norm.counts)
module_trait_cor = cor(module_eigengene,traits,use="p")
module_trait_cor_pvals=corPvalueStudent(module_trait_cor,nSamples)
rownames(traits)=rownames(module_eigengene)
rownames(module_eigengene)



### Visuazlize using heatmap
heatmap_data = merge(module_eigengene,traits,by="row.names")
head(heatmap_data)
heatmap_data = heatmap_data %>%
  column_to_rownames(var="Row.names")

CorLevelPlot(heatmap_data,
             x=names(heatmap_data)[41:44],
             y=names(heatmap_data)[1:40],
             col=c("blue1","skyblue","white","pink","red"))


module_gene_map = as.data.frame(module_colors$colors)
turquoise_genes = module_gene_map %>%
  filter(`module_colors$colors`=="turquoise") %>%
  rownames()

head(turquoise_genes)
gene_symbols_turquoise <- mapIds(
  hgu133plus2.db, 
  keys = turquoise_genes, 
  column = "SYMBOL", 
  keytype = "PROBEID", 
  multiVals = "first"
)

common_de_genes_genesymbol <- mapIds(
  hgu133plus2.db, 
  keys = common_de_genes, 
  column = "SYMBOL", 
  keytype = "PROBEID", 
  multiVals = "first"
)



head(gene_symbols_turquoise)
head(common_de_genes_genesymbol)

# searching for common genes between DE and modules
# Convert to unique gene symbols (to remove duplicates)
unique_gene_symbols_turquoise <- unique(gene_symbols_turquoise)
unique_common_de_genes_genesymbol <- unique(common_de_genes_genesymbol)

# Find common gene symbols
common_genes_de_module <- intersect(unique_gene_symbols_turquoise, unique_common_de_genes_genesymbol)

# Display the common genes
common_genes_de_module

### Highly connected intramdular hub genes
module.membership_measure = cor(module_eigengene,norm.counts,use="p")
module.membership_measure_pvals=corPvalueStudent(module.membership_measure,nSamples)

module.membership_measure_pvals[1:10,1:10]


# Display the common genes
print(common_between_wgcna_and_de)


common_between_wgcna_and_de <- mapIds(hgu133plus2.db, 
                       keys = common_between_wgcna_and_de, 
                       column = "SYMBOL", 
                       keytype = "PROBEID", 
                       multiVals = "first")

result_df <- data.frame(ProbeID = names(common_between_wgcna_and_de), 
                        GeneSymbol = common_between_wgcna_and_de, 
                        stringsAsFactors = FALSE)

# Remove NA values (if any)
result_df <- na.omit(result_df)

# View the resulting gene symbols
print(result_df)

print(result_df$GeneSymbol)

## gene for each module

table(module_colors$colors)

#dendrogram 

plotDendroAndColors(module_colors$dendrograms[[1]]$order, cbind(module_colors$unmergedColors, module_colors$colors),c("unmerged","merged"),dendroLabels = FALSE,addGuide = TRUE,
                    hang=0.03,guideHang = 0.05)


#############checks 

cat("Number of genes in dendrogram:", length(module_colors$dendrograms[[1]]$order), "\n")
cat("Number of genes in unmergedColors:", length(module_colors$unmergedColors), "\n")
cat("Number of genes in colors:", length(module_colors$colors), "\n")

# Module analysis
module_eigengene <- module_colors$MEs
print("Module sizes:")
print(table(module_colors$colors))

# Plot dendrogram with proper color matching
colors_to_plot <- data.frame(
  unmerged = module_colors$unmergedColors[module_colors$dendrograms[[1]]$order],
  merged = module_colors$colors[module_colors$dendrograms[[1]]$order]
)

plotDendroAndColors(
  module_colors$dendrograms[[1]], 
  colors_to_plot,
  c("unmerged", "merged"),
  dendroLabels = FALSE,
  addGuide = TRUE,
  hang = 0.03,
  guideHang = 0.05
)

save.image('D:\\projects\\DGE\\final_workspace.RData')




