library(ggplot2)

# Load GSEA data

gsea_plot <- read_csv("gsea_plot.csv")

# Set variable levels for plotting

recist_levels <-c("PD","SD","PR","CR")
signature_levels <- c("IMM_RES","DC3","DC2","ACTIVATED_DC","IMMATURE_DC","TNFA","DC1","EXHAUSTION","TPEX")

# Create plot

ggplot(gsea_plot, aes(x = factor(RECIST, levels = recist_levels),
                      y = factor(Signature, levels = signature_levels))) + 
geom_point(aes(size = neglog10_FDR, color = NES)) + 
theme_bw() +
theme(panel.grid.major = element_line()) +
scale_color_gradientn(colours = c("blue", "white", "red")) +
scale_x_discrete(position = "bottom") +
labs(x = NULL, y = NULL) +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank()) +
facet_grid(~ Response, scales = "free_x", space = "free_x")
