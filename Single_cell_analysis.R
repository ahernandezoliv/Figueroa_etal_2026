set.seed(1234)
library(Seurat)
library(clustree)
library(dittoSeq)
library(tidyverse)
library(ggpubr)
library(smplot2)

# Loading of raw counts and Seurat object creation from pre-treatment samples

counts_p1 <- read.table("patient1_T0_CD45-counts.tsv")
p1 <- CreateSeuratObject(counts = "counts_p1")
t0_subset <- merge(x = p1, y = c(all other patients))
t0_subset[["RNA"]] <- JoinLayers(t0_subset)[["RNA"]]

# Data processing and integration. This script was applied to the analysis of all cells as well as to the analysis of CD8+ T cells, Myeloid cells, and tumor cell subsets and integration. The number of dimensions  
# used during integration, neighbors, and UMAP analysis was modyfied according to each analysis

t0_subset[["RNA"]] <- split(t0_subset[["RNA"]], f = t0_subset$Patient)
t0_subset <- NormalizeData(t0_subset)
t0_subset <- FindVariableFeatures(t0_subset)
t0_subset <- ScaleData(t0_subset)
t0_subset <- RunPCA(t0_subset, npcs = 100)
ElbowPlot(t0_subset, ndims = 100)
t0_subset <- IntegrateLayers(object = t0_subset, 
                             method = CCAIntegration, 
                             orig.reduction = "pca", 
                             new.reduction = "integrated.cca",
                             dims = 1:75)       
t0_subset[["RNA"]] <- JoinLayers(t0_subset[["RNA"]])
t0_subset <- FindNeighbors(t0_subset,
                           reduction = "integrated.cca", 
                           dims = 1:75)
t0_subset <- FindClusters(t0_subset, 
                          resolution = seq(0.1, 2, by = 0.1),
                          algorithm = 4, random.seed = c(1234))
clustree(t0_subset, prefix = "RNA_snn_res.")
t0_subset <- RunUMAP(t0_subset, dims = 1:75, reduction = "integrated.cca")

# Signature analysis for cluster annotation. Genes were taken from the original paper

tumor_markers <- list(c("SERPINE2","S100A13","MLANA","S100B","MITF"))
naive_markers <- list(c("LTB","CCR7","SELL","IL7R","TCF7","FOS","LEF1","GPR183","KLF2","SLC2A3",
                        "FLT3LG","TPT1","NOSIP","EEF1B2","FXYD5","JUNB","EEF1G","TMEM123","PABPC1","SPOCK2"))
em_like_markers <- list(c("ITM2C","GZMK","GZMA","GNLY","CCL5","GZMM","GZMH","CXCR4","ALOX5AP","ANXA1",
                          "AHNAK","CD52","LGALS3","CD8B","TC2N","C12orf75","ZFP36L2","AOAH","CD8A","PARP8"))
tpex_markers <- list(c("XCL2","CCL4L2","CXCL13","XCL1","TNFRSF9","FABP5","PKM","CRTAM","GAPDH","STMN1","TUBA1B",
                       "TUBB","MIR155HG","CD82","CCL4","TIGIT","PARK7","RAN","ZBED2","TPI1"))
texh_markers <- list(c("VCAM1","HAVCR2","CCL3","LAG3","LINCO1943","TNFRSF9","NKG7","DUSP4","GZMB",
                       "CTLA4","CD27","PHLDA1","CXCR6","CXCL13","PRF1","TIGIT","ADGRG1","ENTPD1","PDCD1","ACP5"))
cx3cr1_markers <- list(c("FGFBP2","CX3CR1","GZMH","FCGR3A","PLEK","GNLY","KLRD1","GZMB","SPON2","NKG7","S1PR5",
                         "FLNA","C1orf21","PRF1","TGFBR3","KLRG1","BIN2","EMP3","S100A4","GZMM"))
hsp_markers <- list(c("HSPA1A","HSPA6","HSPA1B","DNAJB1","HSP90AA1","HSPH1","HSPE1","HSPD1","MTRNR2L8",
                      "HSPB1","HSP90AB1","CACYBP","ZFAND2A","BAG3","DNAJA1","UBC","HSPA8","NR4A1","DUSP1",
                      "DNAJB4"))
isg_markers <- list(c("ISG15","MX1","IFI6","IFIT3","IFIT1","OAS1","RSAD2","IFI44L","IFIT2","MX2","ISG20",
                      "OAS3","STAT1","EPSTI1","IRF7","EIF2AK2","LY6E","IFI35","IFITM1","XAF1"))
pdc_markers <- list(c("GZMB","JCHAIN","PTGDS","IRF7","ITM2C","LILRA4","CLIC3","LTB","TSPAN13",
                      "IRF8","IRF4","PLD4","C12orf75","TCL1A","IL3RA","SPIB","DERL3","CCDC50",
                      "MZB1","PLAC8"))
dc1_markers <- list(c("IDO1","CPVL","CPNE3","DNASE1L3","LGALS2","CLEC9A","TXN","IRF8","CCND1",
                      "CST3","C1orf54","SNX3","S100B","WDFY4","LTB","TAP1","GSTP1","NAPSA","ASB2","NAAA"))
dc2_markers <- list(c("HLA-DQB2","S100B","CD1A","CD207","FCER1A","LTB","TACSTD2","CD1E","NDRG2","CD1C",
                      "FCGBP","ACOT7","GSN","HLA-DQA2","HLA-DQA1","PLEK2","SYNGR2","SERPINF1","C15orf48",
                      "CST3"))
dc3_markers <- list(c("CCL17","CCL22","BIRC3","CCR7","LAMP3","FSCN1","CCL19","TXN","MARCKSL1",
                      "IDO1","EBI3","IL7R","KIF2A","NUB1","CRIP1","LY75","DUSP4","TRAF1","LAD1","KDM2B"))
mac_s100a8_markers <- list(c("S100A8","S100A9","S100A12","THBS1","FCN1","VCAN","EREG","IL1B","CCL20",
                             "LYZ","TIMP1","CXCL3","NLRP3","CXCL2","G0S2","MCEMP1","CD300E","PLAUR","APOBEC3A","ATP2B1-AS1"))
mac_trem2_markers <- list(c("SPP1","PLIN2","FN1","CXCL3","TREM2","CCL3","MT1X","CCL4L2","CCL4","CCL3L1","CCL2","FABP5",
                            "FCGBP","HAMP","MIF","CSTB","LGALS1","CD14","C3","APOE"))
monodc_markers <- list(c("CD1C","FCER1A","HSPA6","CD1E","CLEC10A","HLA-DQA1","DNAJB1","HSPA1A","BAG3","AREG","HSPA1B",
                         "HSPB1","G0S2","PKIB","ZFAND2A","CST7","GSN","CD1A","S100B","HSPH1"))
mono_cd16_markers <- list(c("IFITM2","SMIM25","LST1","CFD","LILRA5","LILRB2","PECAM1","WARS","CFP",
                            "CD52","FCN1","TCF7L2","IRF1","IFITM3","STXBP2","SERPINA1","FCGR3A","CTSS","APOBEC3A","MTSS1"))
mac_isg_markers <- list(c("ISG15","CXCL10","IFIT1","IFIT3","MX1","IFIT2","RSAD2","IFI6","IFITM1","LY6E","APOBEC3A",
                          "TNFSF10","IFITM3","IFI27","IFI44L","ISG20","MX2","OAS1","OAS3","HERC5"))
mac_c1q_markers <- list(c("SELENOP","RNASE1","SLC40A1","FOLR2","PLTP","C1QB","IFI27","C1QA","APOE",
                          "APOC1","NUPR1","C1QC","LGMN","PLD3","DAB2","CTSD","FTL","GPNMB","A2M","F13A1"))
mac_cxcl9_markers <- list(c("CXCL10","CXCL9","GBP1","GBP5","CXCL11","WARS","IDO1","GBP4","CALHM6","VAMP5",
                            "PLAAT4","STAT1","SERPING1","TNFSF10","LAP3","ANKRD22","APOBEC3A","IRF1","PSME2","SLAMF8"))

t0_subset <- AddModuleScore(t0_subset,
                            features = tpex_markers,
                            name = "Tpex")

# Plots

# Dot plot
final_markers <- c("SERPINE2","S100A13","MLANA",
                   "PTPRC","CD3E","CD8A",
                   "CCR7","SELL","IL7R",
                   "ITM2C","GZMK","GZMA",
                   "XCL2","XCL1","CD200",
                   "VCAM1","HAVCR2","CCL3",
                   "FGFBP2","CX3CR1","GZMH",
                   "HSPA1A","HSPA6","HSPA1B",
                   "ISG15","MX1","IFI6",
                   "CD4","CD68","CD86","ITGAM","ITGAX","CD14","FCGR3A",
                   "CLEC9A","XCR1","BATF3",
                   "CCL22","LAMP3","CCL19",
                   "S100A8","S100A9","S100A12",
                   "TREM2","C1QB","C1QA",
                   "CD1C","FCER1A","CD1E")

signatures <- c("Tumor","Naive-like","EM-like","Tpex","Texh",
                "CX3CR1+","HSP","ISG","cDC1",
                "cDC3","Mac-S100A8","Mac-TREM2",
                "Mac-C1Q","MonoDC")

DotPlot(t0_subset, features = c(final_markers,
                                signatures), 
        col.min = 0,dot.min = 0.1,
        group.by = "Final_clusters") + RotatedAxis()+
  geom_vline(xintercept = c(3.5,6.5,27.5,34.5,
                            49.5,63.5), linetype = "solid")+
  geom_hline(yintercept = c(5.5,12.5,13.5), linetype = "solid")+
  ylab("")+xlab("")

# Bar plot. Data was obtained from the metadata ob the seurat object.

ggplot(plot_data, aes(x = Response, y = prop, fill = Final_clusters)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  ylab("Proportion") +
  xlab("Response") +
  labs(fill = "Final Clusters") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  )

# Box plot. Data was obtained from the metadata ob the seurat object.

ggplot(dc1_data, aes(x = Response, y = activ_dc1)) +
  geom_boxplot(outlier.shape = NA,aes(fill = Response), size = 1) +
  stat_compare_means(comparisons = list(c("R", "NR")),
                     method = "t.test",
                     label = "p.signif") +
  scale_fill_manual(values = c("R" = "firebrick",
                               "NR" = "grey"))+
  labs(title = "Activated DC signature in cDC1",
       y = "Activated DC signature score",
       x = "Response") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x  = element_text(size = 14, face = "bold"),
        axis.text.y  = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12))

