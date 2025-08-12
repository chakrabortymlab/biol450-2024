library(dplyr)
library(ggplot2)
 
### R script for enrichment analysis of SVs in marker genes
# generates a null distribution by sampling euchromatic protein coding genes affected by SVs in 10 genomes from this publication (bio450) and DSPR
# Author: Alex Samano, 2025

# SV enrichment 
#KDE is used to inform the rejection sampling process by comparing candidate gene lengths against the empirical distribution (marker genes)

set.seed(111)  # for reproducibility

# load data 
marker_genes <- read.table("marker_mut_genes_SV.txt", header = TRUE, stringsAsFactors = FALSE) 
dspr_sv_genes <- read.table("dmel_proteincoding_euchrom_genes_dsprLSVcount.txt", header = F,col.names = c("gene","length","SV_pres"), stringsAsFactors = FALSE)
bio450_sv_genes <- read.table("dmel_proteincoding_euchrom_genes_biol450LSVcount.txt", header = F,col.names = c("gene","length","SV_pres"), stringsAsFactors = FALSE)

# get number of genes affected by SVs in 10 genomes
sv_genes<-sum(bio450_sv_genes$SV_pres)
sv_genes_percent<-sv_genes/nrow(bio450_sv_genes) * 100

# estimate the density of marker gene lengths
marker_lengths <- marker_genes$length
marker_kde <- density(marker_lengths, kernel = "gaussian", n = 2048) #gaussian or normal

# function to interpolate density at any length
interpolate_density <- approxfun(marker_kde$x, marker_kde$y, rule = 2)

# get observed count of SV affected marker genes
observed_sv_count <- sum(marker_genes$SV_pres)

# run monte carlo simulation 100K times
num_simulations <- 100000

### dspr genome SVs as null
dspr_sv_counts <- numeric(num_simulations)
marker_n <- nrow(marker_genes)

for (i in 1:num_simulations) {
  sampled_genes <- data.frame()
  attempts <- 0
  
  while (nrow(sampled_genes) < marker_n && attempts < 10000) {
    candidate <- dspr_sv_genes[sample(nrow(dspr_sv_genes), 1), ]
    density_at_length <- interpolate_density(candidate$length)
    accept_prob <- density_at_length / max(marker_kde$y)
    
    if (runif(1) < accept_prob && !(candidate$gene %in% sampled_genes$gene)) {
      sampled_genes <- rbind(sampled_genes, candidate)
    }
    attempts <- attempts + 1
  }
  
  dspr_sv_counts[i] <- sum(sampled_genes$SV_pres)
}


### bio450 genomes as null
bio450_sv_counts <- numeric(num_simulations)
marker_n <- nrow(marker_genes)

for (i in 1:num_simulations) {
  sampled_genes <- data.frame()
  attempts <- 0
  
  while (nrow(sampled_genes) < marker_n && attempts < 10000) {
    candidate <- bio450_sv_genes[sample(nrow(bio450_sv_genes), 1), ]
    density_at_length <- interpolate_density(candidate$length)
    accept_prob <- density_at_length / max(marker_kde$y)
    
    if (runif(1) < accept_prob && !(candidate$gene %in% sampled_genes$gene)) {
      sampled_genes <- rbind(sampled_genes, candidate)
    }
    attempts <- attempts + 1
  }
  
  bio450_sv_counts[i] <- sum(sampled_genes$SV_pres)
}


# stats
#calculate mean and median
bio450_median<-median(bio450_sv_counts) # 12
dspr_median<-median(dspr_sv_counts) # 18

# calculate enrichment percent
bio450_enrichment_percent <- ((observed_sv_count - bio450_median) / bio450_median) * 100
dspr_enrichment_percent <- ((observed_sv_count - dspr_median) / dspr_median) * 100

# p value calculation
bio450_p_value <- (sum(bio450_sv_counts >= observed_sv_count) + 1) / (num_simulations + 1)
dspr_p_value <- (sum(dspr_sv_counts >= observed_sv_count) + 1) / (num_simulations + 1)

