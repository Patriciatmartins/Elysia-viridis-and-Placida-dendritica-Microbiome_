##############################################
#
# Title: R Script for Microbiome characterization of the sea slugs Elysia viridis and 
#Placida dendritica: insights into potential roles in kleptoplasty
#
# Description: This script reproduces the analyses and figures presented in the manuscript,
# including Rarefaction (Fig. 1), Venn Diagram (Fig. 2), Alpha Diversity (Fig. 3),
# PCoA & PERMANOVA (Fig. 4), Comparative Taxonomic Composition (Fig. 5),
# and Core Microbiome characterization.
#
## Data availability:
# The input data files (zOTU table and sample metadata) are available
# in the public GitHub repository associated with this study:
# https://github.com/Patriciamartins/Elysia-viridis-and-Placida-dendritica-Microbiome
#
# Files required for running this script:
# - zotus.OTU.table_new.txt  → zOTU abundance and taxonomy table
# - Labels.txt               → sample metadata file
#
#
# Author: Patricia Martins
# Date: March 2025
# R version: 4.2.3
##############################################

# --- Load required libraries ---
library(vegan)
library(VennDiagram)
library(gplots)
library(ggplot2)
library(ggpubr)
library(car)
library(svglite)
library(gridExtra)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(tidyverse)

# --- Set working directory ---
# Before running the script, set your working directory to the folder
# containing the input data files.
# Example (uncomment and modify if needed):
# setwd("path/to/your/folder")

# If you are using RStudio, the following line will automatically set the working 
# directory to the folder where this script is located.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# --- Load zOTU table (Table S3 from the supplementary material) ---
patricia.df <- read.table("zotus.OTU.table_new.txt", header = TRUE, as.is = TRUE, row.names = 1)

# --- Define biological replicates (samples) ---
patricia.mat <- patricia.df[, c("Ev01","Ev02","Ev03","Ev04","Ev05","Ev06","Ev07","Ev08","Ev09","Ev10",
                                "Pd01","Pd02","Pd03","Pd04","Pd05","Pd06","Pd07","Pd08","Pd09","Pd10")]


##################################################
# SECTION 1 — Rarefaction Curves (Fig. 1)
##################################################


otu_table <- t(patricia.mat)
sample_groups <- factor(gsub("[0-9]+", "", rownames(otu_table)))
sample_colors <- ifelse(sample_groups == "Ev", "olivedrab", "orange3")

pdf("Figure 1.pdf", width = 8, height = 6)
rare_data <- rarecurve(otu_table, step = 20, label = FALSE, col = sample_colors, lwd = 2,
                       xlab = "Number of Reads", ylab = "Number of ZOTUs",
                       main = "Rarefaction Curves")
for (i in seq_along(rare_data)) {
  text(max(rare_data[[i]]), rare_data[[i]][length(rare_data[[i]])],
       labels = rownames(otu_table)[i], pos = 4, cex = 0.7, col = sample_colors[i])
}
legend("bottomright", legend = c("E. viridis", "P. dendritica"),
       col = c("olivedrab", "orange3"), lty = 1, lwd = 2, cex = 0.8)
dev.off()


##################################################
# SECTION 2 — Venn Diagram (Fig. 2)
##################################################


patricia.mat_Ev <- patricia.mat[, grep("Ev", colnames(patricia.mat))]
patricia.mat_Pd <- patricia.mat[, grep("Pd", colnames(patricia.mat))]

E.viridis <- names(rowSums(patricia.mat_Ev)[rowSums(patricia.mat_Ev) >= 1])
P.dendritica <- names(rowSums(patricia.mat_Pd)[rowSums(patricia.mat_Pd) >= 1])

svglite("Figure 2.svg", width = 5, height = 5)
draw.pairwise.venn(
  area1 = length(E.viridis),
  area2 = length(P.dendritica),
  cross.area = length(intersect(E.viridis, P.dendritica)),
  category = c(expression(italic("E. viridis")), expression(italic("P. dendritica"))),
  fill = c("olivedrab", "orange3"), alpha = 0.5,
  cex = 1.5, cat.cex = 1.5, cat.pos = c(-20, 20), cat.dist = c(0.05, 0.05)
)
dev.off()


##################################################
# SECTION 3 — Alpha Diversity (Fig. 3)
##################################################


patricia.mat.red2<-patricia.df[,c("Ev01","Ev02","Ev03","Ev04","Ev05", "Ev06", "Ev07","Ev08","Ev09","Ev10","Pd01","Pd02","Pd03","Pd04","Pd05","Pd06","Pd07","Pd08", "Pd09","Pd10")]

species.abundance.red2<-apply(patricia.mat.red2,1,sum)

sample.abundance.red2<-apply(patricia.mat.red2,2,sum)


labels<-read.table("Labels.txt", header=TRUE,as.is=TRUE,row.names=1)

patricia.mean1<-t(patricia.mat.red2)

patricia.mean<-apply(patricia.mean1,2, function(x){tapply(x, labels$Sp, mean)})


patricia.mean <- apply(patricia.mean1, 2, function(x){tapply(x, labels$Sp, mean)})

patricia.sd <- apply(patricia.mean1, 2, function(x){tapply(x, labels$Sp, sd)})


patricia.stats <- data.frame(
  Mean = as.vector(patricia.mean),
  SD = as.vector(patricia.sd)
)
#rownames(patricia.stats) <- rownames(patricia.mean)

shannons.diversity <- diversity(patricia.mean, index = "shannon")
shannons.diversity_ <- diversity(t(patricia.mat.red2), index = "shannon")

Species.richness <- specnumber(patricia.mean)
Species.richness_ <- specnumber(t(patricia.mat.red2))

Pielous.evenness <- shannons.diversity / log(Species.richness)
Pielous.evenness_ <- shannons.diversity_ / log(Species.richness_)

chao1 <- estimateR(t(patricia.mat.red2))["S.chao1", ]

shannons.mean <- tapply(shannons.diversity_, labels$Sp, mean)
shannons.sd <- tapply(shannons.diversity_, labels$Sp, sd)
shannons.stats <- data.frame(Mean = shannons.mean, SD = shannons.sd)

species.mean <- tapply(Species.richness_, labels$Sp, mean)
species.sd <- tapply(Species.richness_, labels$Sp, sd)
species.stats <- data.frame(Mean = species.mean, SD = species.sd)

pielou.mean <- tapply(Pielous.evenness_, labels$Sp, mean)
pielou.sd <- tapply(Pielous.evenness_, labels$Sp, sd)
pielou.stats <- data.frame(Mean = pielou.mean, SD = pielou.sd)

chao1.mean <- tapply(chao1, labels$Sp, mean)
chao1.sd <- tapply(chao1, labels$Sp, sd)
chao1.stats <- data.frame(Mean = chao1.mean, SD = chao1.sd)

Species.richness<-specnumber(patricia.mean)
Species.richness_ <- data.frame(Species.richness)

Species.richness_all<-specnumber(t(patricia.mat.red2))

#Test for significance

diversity_data <- data.frame(
  Group = labels$Sp,  
  Chao1 = chao1,  
  Shannon = shannons.diversity_,  
  Observed_OTU = Species.richness_all,  
  Pielou = Pielous.evenness_  
)

# Test for normality
shapiro.test(diversity_data$Observed_OTU[diversity_data$Group == "Ev"])
shapiro.test(diversity_data$Observed_OTU[diversity_data$Group == "Pd"])

# Test for homogeneity of variances
leveneTest(Observed_OTU ~ Group, data = diversity_data)
t.test(Observed_OTU ~ Group, data = diversity_data, var.equal = TRUE)

shapiro.test(diversity_data$Chao1[diversity_data$Group == "Ev"])
shapiro.test(diversity_data$Chao1[diversity_data$Group == "Pd"])

leveneTest(Chao1 ~ Group, data = diversity_data)
t.test(Chao1 ~ Group, data = diversity_data, var.equal = TRUE)

shapiro.test(diversity_data$Shannon[diversity_data$Group == "Ev"])
shapiro.test(diversity_data$Shannon[diversity_data$Group == "Pd"])

leveneTest(Shannon ~ Group, data = diversity_data)
t.test(Shannon ~ Group, data = diversity_data, var.equal = TRUE)

shapiro.test(diversity_data$Pielou[diversity_data$Group == "Ev"])
shapiro.test(diversity_data$Pielou[diversity_data$Group == "Pd"])

leveneTest(Pielou ~ Group, data = diversity_data)
t.test(Pielou ~ Group, data = diversity_data, var.equal = TRUE)


#Barplot

svglite(file="Figure 3.svg")

opar <- par(no.readonly=TRUE)
par(mfrow = c(1,4), omi=c(3,0.4,0,0), plt=c(0.23,1,0,0.68), cex.main=1, xpd=TRUE)

# Observed OTU
bar1 <- barplot2(height=as.matrix(Species.richness_[,1]),las=0,names.arg=c("Ev","Pd"), main= "Observed OTU", beside=TRUE, space=c(0.25,0.25),  ci.l = Species.richness_[,1], ci.u = Species.richness_[,1] + Species.richness_[,2], col=c("olivedrab", "orange3"),ylim=c(0,200))

# Chao 1
bar2 <- barplot2(height=as.matrix(chao1.stats[,1]),
                 las=0, names.arg=c("Ev", "Pd"), main="Chao 1", beside=TRUE, 
                 space=c(0.25, 0.25), plot.ci=TRUE, ci.l=chao1.stats[,1], 
                 ci.u=chao1.stats[,1] + chao1.stats[,2], 
                 col=c("olivedrab", "orange3"), ylim=c(0,100))
text(x=bar2, y=chao1.stats[,1] + 5, 
     labels=c("t = -4.48, p < 0.001", ""), cex=1, col="black")  # Adjust values

# Shannon index
bar3 <- barplot2(height=as.matrix(shannons.stats[,1]),
                 las=0, names.arg=c("Ev", "Pd"), main="Shannon index", beside=TRUE, 
                 space=c(0.25, 0.25), plot.ci=TRUE, ci.l=shannons.stats[,1], 
                 ci.u=shannons.stats[,1] + shannons.stats[,2], 
                 col=c("olivedrab", "orange3"), ylim=c(0,4))
text(x=bar3, y=shannons.stats[,1] + 0.2, 
     labels=c("t = -7.69, p < 0.001", ""), cex=1, col="black")  # Adjust values

# Pielou's evenness
bar4 <- barplot2(height=as.matrix(pielou.stats[,1]),
                 las=0, names.arg=c("Ev", "Pd"), main="Pielou's evenness", beside=TRUE, 
                 space=c(0.25, 0.25), plot.ci=TRUE, ci.l=pielou.stats[,1], 
                 ci.u=pielou.stats[,1] + pielou.stats[,2], 
                 col=c("olivedrab", "orange3"), ylim=c(0,1))
text(x=bar4, y=pielou.stats[,1] + 0.05, 
     labels=c("t = -6.92, p < 0.001", ""), cex=1, col="black")  # Adjust values

dev.off()


##################################################
# SECTION 4 — PCoA Ordination & PERMANOVA (Fig. 4)
##################################################

patricia.bray.dist <- vegdist(t(log(patricia.mat + 1)), method = "bray")
patricia.pco <- cmdscale(patricia.bray.dist, eig = TRUE)
patricia.wa1 <- wascores(patricia.pco$points, t(patricia.mat))
species.abundance <- apply(patricia.mat, 1, sum)

# --- PCoA ordination plot ---
svglite("Figure 4.svg", width = 6, height = 6)
translucentGreyBorder <- adjustcolor("grey80", alpha.f = 0.5)
translucentGreyFill   <- adjustcolor("grey80", alpha.f = 0.5)
sample_colors <- c(rep("olivedrab", 10), rep("orange3", 10))
sample_symbols <- rep(22, 20)
x_range <- range(c(patricia.wa1[, 1], patricia.pco$points[, 1]))
y_range <- range(c(patricia.wa1[, 2], patricia.pco$points[, 2]))
plot(x_range, y_range, type = "n", xlab = "Axis 1", ylab = "Axis 2", asp = 1)

groups <- list(
  which(species.abundance < 1000),
  which(species.abundance >= 1000 & species.abundance < 2000),
  which(species.abundance >= 2000 & species.abundance < 3000),
  which(species.abundance >= 3000 & species.abundance < 10000),
  which(species.abundance >= 10000)
)
sizes <- c(0.3, 0.7, 1.0, 1.5, 2.0)
for (i in seq_along(groups)) {
  points(patricia.wa1[groups[[i]], 1], patricia.wa1[groups[[i]], 2],
         pch = 21, cex = sizes[i] * 0.7, col = translucentGreyBorder, bg = translucentGreyFill)
}
points(patricia.pco$points[, 1], patricia.pco$points[, 2],
       pch = sample_symbols, col = sample_colors, bg = sample_colors)
legend("topright", legend = c("E. viridis", "P. dendritica"),
       pch = 22, pt.bg = c("olivedrab", "orange3"), col = c("olivedrab", "orange3"), cex = 0.8)
legend("bottomright", legend = c("< 1000", "1000–2000", "2000–3000", "3000–10000", "> 10000"),
       pch = 21, pt.cex = sizes, pt.bg = rep(translucentGreyFill, 5), col = rep(translucentGreyBorder, 5))
dev.off()

# --- PERMANOVA ---
patricia.mat.t <- t(patricia.mat)
dist_matrix <- vegdist(patricia.mat.t, method = "bray")
groups <- factor(labels$Sp)
permanova_result <- adonis2(dist_matrix ~ groups)
print(permanova_result)



##################################################
# SECTION 5 — Comparative Taxonomic Composition (Fig. 5)
##################################################



patricia.df.red <- patricia.df

patricia.mat2 <- patricia.df[, c("Ev01","Ev02","Ev03","Ev04","Ev05", 
                                 "Ev06","Ev07","Ev08","Ev09","Ev10",
                                 "Pd01","Pd02","Pd03","Pd04","Pd05",
                                 "Pd06","Pd07","Pd08","Pd09","Pd10")]

normalize_relative_abundance <- function(mat) {
  sweep(mat, 2, colSums(mat), FUN = "/") * 100
}

process_taxa_data_with_threshold <- function(mat, level_name, threshold = 1) {
  mat_rel_abund <- normalize_relative_abundance(mat)
  df <- as.data.frame(mat_rel_abund)
  df$Taxon <- rownames(df)
  
  df_long <- df %>%
    pivot_longer(cols = -Taxon, names_to = "Sample", values_to = "RelativeAbundance") %>%
    mutate(Sp = substr(Sample, 1, 2))  
  
  df_summary <- df_long %>%
    group_by(Taxon, Sp) %>%
    summarise(Mean = mean(RelativeAbundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(Level = level_name, Percentage = Mean)
  
  df_summary <- df_summary %>%
    mutate(is_low_abundance = Percentage < threshold) %>%
    mutate(Taxon = ifelse(is_low_abundance, "Other", Taxon)) %>%
    group_by(Sp, Taxon, Level) %>%
    summarise(Percentage = sum(Percentage), .groups = "drop")
  
  return(df_summary)
}

patricia.mat.Phylum <- apply(patricia.mat2, 2, function(x) tapply(x, patricia.df.red$Phylum, sum))
patricia.mat.Class  <- apply(patricia.mat2, 2, function(x) tapply(x, patricia.df.red$Class, sum))
patricia.mat.Order  <- apply(patricia.mat2, 2, function(x) tapply(x, patricia.df.red$Order, sum))

threshold_value <- 1
phyla_data <- process_taxa_data_with_threshold(patricia.mat.Phylum, "Phylum", threshold = threshold_value)
class_data <- process_taxa_data_with_threshold(patricia.mat.Class, "Class", threshold = threshold_value)
order_data <- process_taxa_data_with_threshold(patricia.mat.Order, "Order", threshold = threshold_value)
all_taxa_data <- bind_rows(phyla_data, class_data, order_data)

phylum_colors <- colorRampPalette(brewer.pal(8, "Set1"))
class_colors  <- colorRampPalette(brewer.pal(8, "Set2"))
order_colors  <- colorRampPalette(brewer.pal(8, "Set3"))

# Create stacked bar plots
create_taxa_plot_stacked <- function(taxa_level, title, color_palette) {
  taxa_data <- all_taxa_data %>% filter(Level == taxa_level, Sp %in% c("Ev", "Pd"))
  
  ggplot(taxa_data, aes(x = Sp, y = Percentage, fill = Taxon)) +
    geom_bar(stat = "identity", position = "fill", width = 0.6) +
    scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      labels = function(x) paste0(x * 100, "%")  # Adds "%" to y-axis tick labels
    ) +
    labs(title = title, y = "Relative Abundance", x = NULL) +
    theme_minimal() +
    scale_fill_manual(values = color_palette(length(unique(taxa_data$Taxon)))) +
    theme(
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      plot.title = element_text(size = 12, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 10)
    )
}

plot_phylum <- create_taxa_plot_stacked("Phylum", "A) Phyla", phylum_colors)
plot_class  <- create_taxa_plot_stacked("Class", "B) Class", class_colors)
plot_order  <- create_taxa_plot_stacked("Order", "C) Order", order_colors)

svglite("Figure 5.svg", width = 10, height = 12)
grid.arrange(plot_phylum, plot_class, plot_order, ncol = 1)
dev.off()


  
  ##################################################
  # SECTION 6 — Core Microbiome Analysis
  ##################################################
  
  # Identify zOTUs detected in all replicates of each species
  Ev.samples <- paste0("Ev", sprintf("%02d", 1:10))
  Pd.samples <- paste0("Pd", sprintf("%02d", 1:10))
  
  patricia.mat.Ev <- patricia.df[, Ev.samples]
  patricia.mat.Pd <- patricia.df[, Pd.samples]
  
  present_in_all.Ev <- apply(patricia.mat.Ev, 1, function(x) all(x > 0))
  present_in_all.Pd <- apply(patricia.mat.Pd, 1, function(x) all(x > 0))
  
  core_zOTUs_Ev <- rownames(patricia.mat.Ev)[present_in_all.Ev]
  core_zOTUs_Pd <- rownames(patricia.mat.Pd)[present_in_all.Pd]
  
  Ev_core_perc <- prop.table(as.matrix(Ev_core[, Ev.samples]), 2) * 100
  Pd_core_perc <- prop.table(as.matrix(Pd_core[, Pd.samples]), 2) * 100
  
  write.table(Ev_core_perc, "E_viridis_core_OTUs_relative_abundance.txt", sep = "\t", quote = FALSE)
  write.table(Pd_core_perc, "P_dendritica_core_OTUs_relative_abundance.txt", sep = "\t", quote = FALSE)
  
  Ev_phylum_mean <- tapply(rowMeans(Ev_core_perc), Ev_core$Phylum, mean)
  Pd_phylum_mean <- tapply(rowMeans(Pd_core_perc), Pd_core$Phylum, mean)
  
  barplot(Ev_phylum_mean, col = "olivedrab", main = "Core Phyla - E. viridis",
          ylab = "Mean Relative Abundance (%)", las = 2)
  barplot(Pd_phylum_mean, col = "orange3", main = "Core Phyla - P. dendritica",
          ylab = "Mean Relative Abundance (%)", las = 2)
  
  