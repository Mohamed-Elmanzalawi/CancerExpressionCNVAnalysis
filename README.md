# Computational Analysis of Lung and Kidney Cancer Gene Expressions and CNV Relationships

## Project Overview

This repository documents a collaborative project undertaken by Mohamed El-Manzalawi, Amany Awad, Rania Alanany, and Ahmed Elghamry at Nile University, Egypt. The project, conducted as part of the CIT660: Statistical Analysis and Visualization course, is titled "Computational Analysis of Lung and Kidney Cancer Gene Expressions and their CNV Relationship." The team was under the supervision of Prof. Ibrahim Mohamed Youssef, PhD, from the Faculty of Engineering, Cairo University.

## Project Description

Lung and kidney cancers are the leading causes of cancer-related deaths worldwide. This project aims to analyze gene expressions in Lung Squamous Cell Carcinoma (LUSC) and Kidney Renal Clear Cell Carcinoma (KIRC). The primary objectives include identifying differentially expressed genes (DEGs), exploring copy number variations (CNVs), and performing Gene Set Enrichment Analysis (GSEA) to understand the biological implications.

## Key files Structure

- **Code:** R scripts for data filtration, hypothesis testing, GSEA, CNV filtration, regression, and David enrichment analysis.
  
- **Data:** Raw unprocessed data files, including gene expression data, and CNV data.

- **Report:** Comprehensive report detailing project background, objectives, methodology, results, and discussions.
  
## 1. Introduction

### Lung Cancer: A Global Challenge
Lung cancer, the leading cause of cancer-related death worldwide, poses a significant public health challenge. We focus on Lung Squamous Cell Carcinoma (LUSC) and delve into the genetic factors influencing its development and progression.

### Kidney Cancer: An Enigmatic Foe
Kidney Renal Clear Cell Carcinoma (KIRC), the eighth most common cancer, presents unique challenges due to its resistance to conventional treatments. Our study explores the genomic alterations and molecular mechanisms associated with KIRC.

## 2. Data Collection
Our primary source of genetic data stems from The Cancer Genome Atlas (TCGA) data portal, a comprehensive repository that provides a wealth of information on various cancer types. Specifically, we focused on Lung Squamous Cell Carcinoma (LUSC) and Kidney Renal Clear Cell Carcinoma (KIRC), acknowledging their prominence in cancer-related mortality.

## 3. Data Filtering
The filtration process involved eliminating genes with more than 50% zero values in both cancer and healthy groups for each cancer type (KIRC and LUSC). This not only streamlined subsequent analyses but also ensured the statistical robustness of our findings. The table below summarizes the gene and sample count changes resulting from this filtration process.

![](/Pictures/genes&samples_num.png)

## 4. Methods

### 4.1. Hypothesis Testing
We employ hypothesis testing, fold change, and volcano plots to identify DEGs in both cancer types. Data filtration ensures robust results, and statistical methods guide our analysis.

### 4.2. GSEA
Gene Set Enrichment Analysis (GSEA) adds depth to our study, unraveling the intricate pathways associated with DEGs. The results provide valuable insights into the molecular signatures of lung and kidney cancers.

### 4.3. Regression
Regression analysis establishes connections between gene expressions and CNVs. We utilize the glmnet package for efficient analysis, contributing to a comprehensive understanding of genetic interactions.

### 4.4. Enrichment Analysis
We employ DAVID network analysis to explore enriched Gene Ontology terms, OMIM diseases, and KEGG pathways, enhancing our understanding of the functional significance of DEGs.

## 5. Literature Review

A review of existing literature validates our findings. Genes such as FBXO45, FGF11, and SERPINH1 emerge as potential biomarkers in lung and kidney cancers, aligning with previous research.

## Instructions for Reproduction

To reproduce the analysis:
1. Clone the repository to your local machine:
   ```bash
   git clone https://github.com/your-username/lung-kidney-cancer-analysis.git
2. Open RStudio or any R environment.
3. Execute the R scripts in the Code directory in the provided order to perform data analysis But be sure to change the file paths according to your data location.
4. Refer to the report in the Report folder for detailed results and discussions.

## Contributors

- Mohamed Elmanzalawi [Linkedin](https://www.linkedin.com/in/mohamed-elmanzalawi/)
- Ahmed Nabil Elghamry [Linkedin](https://www.linkedin.com/in/ahmed-elghamry-7b22829a/)
- Amany Awad
- Rania Alanany

## Acknowledgments
We express our sincere gratitude to Prof.Ibrahim Mohamed Youssef [Linkedin](https://https://www.linkedin.com/in/ibrahim-youssef-65262a145/) for his invaluable supervision and guidance throughout the project.

## Keywords
Lung Cancer, Kidney Cancer, Gene Expression, GSEA, CNV, Computational Analysis.


