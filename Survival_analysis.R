set.seed(1234)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(survival)
library(survminer)

# Loading and gene symbol annotation of TCGA-SKCM data. Genes with more than one ensembl id were aggregated according to the mean expression.

skcm_counts <- read.delim("TCGA-SKCM.star_fpkm-uq.tsv.gz",
                        row.names = 1,
                        header = T)

ensembl_ids <- sapply(strsplit(rownames(skcm_counts), 
                               "\\."), `[`, 1)

gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

expr_anno <- data.frame(Symbol = gene_symbols, 
                        skcm_counts, 
                        check.names = FALSE)

expr_anno <- expr_anno[!is.na(expr_anno$Symbol), ]

expr_anno_agg <- expr_anno %>%
  group_by(Symbol) %>%
  summarise(across(where(is.numeric), mean))

expr_anno_agg <- as.data.frame(expr_anno_agg)

rownames(expr_anno_agg) <- expr_anno_agg$Symbol

expr_anno_agg$Symbol <- NULL

# Signature calculation

signatures <- read.csv("signatures_melanoma.csv",
                       header = T)

signature_means <- sapply(colnames(signatures), function(sig) {
  genes <- signatures[[sig]]
  genes_present <- genes[genes %in% rownames(expr_anno_agg)]  
  if (length(genes_present) > 0) {
    mean_values <- colMeans(expr_anno_agg[genes_present, , drop = FALSE], 
                            na.rm = TRUE)  
    mean_values[is.nan(mean_values)] <- NA  
    return(mean_values)
  } else {
    return(rep(NA, ncol(expr_anno_agg)))
  }
})

signature_means <- as.data.frame(signature_means)

# Merge with clinical data

survival_skcm <- read.delim("TCGA-SKCM.survival.tsv.gz", 
                            row.names = 1, 
                            header = T)

survival_skcm <- survival_skcm[,-3]
survival_skcm$Patient <- rownames(survival_skcm)

signature_means$Patient <- rownames(signature_means)

surv_sign <- left_join(survival_skcm, 
                       signature_means, 
                       by = "Patient")

# Survival analysis

variables_skcm <- colnames(surv_sign)[4:37]

data_os <- surv_cutpoint(skcm_data,
                         time = "OS_time",
                         event = "OS",
                         variables = variables_skcm,
                         minprop = 0.25,
                         progressbar = TRUE)     

DataCat_os<- surv_categorize(data_os, labels = c("Low","High"))

theme_survminer_cent <- function() {
  theme_survminer() %+replace%
    theme(plot.title = element_text(size = 14,hjust = 0.5),
          plot.subtitle = element_text(size = 13,hjust = 0.5),
          legend.key.size = unit(0.6, "lines"),
          legend.margin = margin(t = -7, unit = "points")
    )
}

kmPlot_os <- ggsurvplot(surv_fit(Surv(OS_time, OS)~Signature, data = DataCat_os),
                        data = DataCat_os,
                        risk.table = F,
                        conf.int = F,
                        pval = TRUE,
                        pval.method = TRUE,
                        xscale = "d_y",
                        surv.scale = "percent",
                        xlab = "Time to event (years)",
                        ylab = "Survival probability (%)",
                        break.time.by = 1826,
                        legend = "bottom",
                        legend.title = "",
                        legend.labs = c("High","Low"),
                        title = "Signature",
                        censor.size = 2,
                        size = 0.5,
                        palette = c("red", "blue"),
                        ggtheme = theme_survminer_cent())
                          
kmPlot_os
