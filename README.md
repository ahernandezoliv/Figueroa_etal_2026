# Figueroa_etal_2026
This repository contain the script used for the analysis and plotting of the data from the paper: "Lymphodepleting preconditioning impairs host antitumor immunity induced by adoptive cell therapy in mouse models"

This repository is divided in three different R scripts:

1. GSEA analysis. Includes the script for generating the dotplot in Figure 7a. No other script was used for this plot science normalized data was download directly from GEO and GSEA (Broad Institute) is graphic tool that uses no code.
2. Single-cell analysis. Includes the script for processing and integration of single-cell data with corresponding plots available at Figures 7b-7d and Supplementary Figure 13.
3. Survival analysis. Includes the script for calculating the mean expression of each gene signature and subsequently survival analysis and Kaplan-Meier curves present in Figure 7e.

Data used in this study:

1. Bulk and single-cell data of ACT-treated patients were obtained from Gene Expression Omnibus with accession codes GSE100797 and GSE222448, respectively.
2. TCGA-SKCM gene expression and clinical data were downloaded from Xena Browser (https://xenabrowser.net/datapages/).
