### Idiopathic Interstitial Pneumonia (IIP) Gene Expression Analysis

This repository contains the code and workflow used for analyzing gene expression data from the GSE32537 dataset. The goal of the project is to explore the molecular phenotyping of Idiopathic Interstitial Pneumonias (IIPs), identify differentially expressed genes (DEGs), perform Gene Set Enrichment Analysis (GSEA), and conduct Weighted Gene Co-expression Network Analysis (WGCNA) on the dataset. This is a preliminary workflow, and further biological interpretation of the identified modules is ongoing.
Dataset Information

    Dataset Name: GSE32537
    Status: Public (Available since June 21, 2013)
    Organism: Homo sapiens
    Experiment Type: Expression profiling by array
    Summary:
        The dataset consists of transcriptional and miRNA profiles from lung tissue of 167 subjects with idiopathic interstitial pneumonias (IIP) and 50 non-diseased controls.
        Goal: To investigate gene expression patterns across different clinical subtypes of IIP to identify molecular signatures associated with disease severity and prognosis.
        Findings: A strong molecular signature associated with cilium genes was identified, linking these genes to better survival outcomes in IIP patients.

### Analysis Workflow

    Differential Gene Expression (DGE) Analysis:
        Performed DGE analysis to identify common DE genes across various severity levels of IIP.

    Gene Set Enrichment Analysis (GSEA):
        Conducted GSEA to identify enriched biological pathways among the DE genes.

    Weighted Gene Co-expression Network Analysis (WGCNA):
        Used WGCNA to identify gene modules and associated hub genes.
        Focused on identifying common genes between the DEGs and selected WGCNA module genes.

    Preliminary Results:
        The workflow has identified common genes between DEGs and certain gene modules; however, biological relevance and interpretation of the identified modules are still under investigation.

### Dependencies

The following R libraries are required to run the analysis. The installation commands for each library are provided below.

# Load the required libraries
library(impute)
library(WGCNA)
library(biomaRt)
install.packages("doParallel")
install.packages("foreach")
library(doParallel)
library(foreach)
install.packages("BiocManager")
BiocManager::install("beadarray")
library(CorLevelPlot)

# Load remaining required libraries
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
BiocManager::install("clusterProfiler", force = TRUE)
BiocManager::install("arrayQualityMetrics")
BiocManager::install("VennDiagram")
library(arrayQualityMetrics)
install.packages("splitstackshape")
library("splitstackshape")
library("tidyr")
library(limma)
library(org.Hs.eg.db)
library(dplyr)
library(VennDiagram)

if (!requireNamespace("hugene10sttranscriptcluster.db", quietly = TRUE)) {
  BiocManager::install("hugene10sttranscriptcluster.db")
}
library(hugene10sttranscriptcluster.db)

## Repository Structure

The repository contains only the `IPF_severity_DGE_and_WGCNA_flow.R` script, which performs the differential gene expression (DGE) analysis and Weighted Gene Co-expression Network Analysis (WGCNA).  
- **/modified_target.csv**: A CSV file containing sample IDs and conditions used as input for the analysis.  


How to Run the Code

1. Clone this repository to your local machine:

    ```bash
    git clone https://github.com/yourusername/repository-name.git
    ```

2. Ensure you have the necessary R packages installed. Use the commands provided in the script or the dependency installation section.

3. Follow these steps in R:
    - Load the `modified_target.csv` file and `GSE32537_family.soft` dataset.
    - Run the `IPF_severity_DGE_and_WGCNA_flow.R` script to:
        - Perform Differential Gene Expression (DGE) analysis.
        - Conduct Weighted Gene Co-expression Network Analysis (WGCNA).




