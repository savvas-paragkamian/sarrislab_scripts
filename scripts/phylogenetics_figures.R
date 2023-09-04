#!/usr/bin/env Rscript

###############################################################################
# script name: phygenetics_figures.R
# developed by: Savvas Paragkamian
# framework: SarrisLab
###############################################################################
# GOAL:
# Aim of this script is to analyse and  visualise the de novo trees from the
# gtdb-tk pipeline.
###############################################################################
# usage:./phygenetics_figures.R
###############################################################################

library(tidyverse)
library(ape)
library(ggtree)
library(tidytree)

# read the taxonomy
taxonomy <- read_delim( "../543_de_novo/infer/gtdbtk.taxonomy_long.tsv",
                       delim="\t",
                       col_names=T)
# transform taxonomy to wide format
taxonomy_w <- taxonomy %>% group_by(gtdb_id, classification) %>%
    summarise(taxa=paste(taxonname, collapse = "; "),
              .groups = "drop") %>%  
    pivot_wider(names_from=classification, values_from=taxa)

# load tree
tree_543 <- ape::read.tree(file = "../543_de_novo/infer/gtdbtk.bac120.decorated.tree")
tree_543_taxonomy <- as_tibble(tree_543) %>% left_join(taxonomy_w, by=c("label"="gtdb_id"))

# the node 82141 is the family
claude_f__DSM_18226 <- groupClade(as_tibble(tree_543), 82141) %>%
    filter(group==1) %>% 
    left_join(taxonomy_w,
              by=c("label"="gtdb_id")) %>% 
    as.treedata()

tree_543_long <- as_tibble(tree_543)
as.treedata(tree_543_long)

# visualisation

tree_plot <- ggtree(data = claude_f__DSM_18226,
               aes(color=g)) + 
     geom_tiplab(size=2, aes(label= s, color=g)) + 
     theme_tree2(legend.position = c(.1, .88))
    
ggsave(plot=tree_plot,
       "../tree_plot1.png",
       device="png",
       height = 30,
       width=50,
       units="cm")
