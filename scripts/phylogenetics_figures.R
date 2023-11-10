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
taxonomy <- read_delim( "../Assemblies_de_novo/infer/gtdbtk.taxonomy_long.tsv",
                       delim="\t",
                       col_names=T)
# transform taxonomy to wide format
taxonomy_w <- taxonomy %>% group_by(gtdb_id, classification) %>%
    summarise(taxa=paste(taxonname, collapse = "; "),
              .groups = "drop") %>%  
    pivot_wider(names_from=classification, values_from=taxa)

# load tree
tree <- ape::read.tree(file = "../Assemblies_de_novo/infer/gtdbtk.bac120.decorated.tree")
tree_taxonomy <- as_tibble(tree) %>% left_join(taxonomy_w, by=c("label"="gtdb_id"))

# filter the tree

as_tibble(tree) %>% filter(grep("Bacilli",label))

# it doesn't have a node for Bacilli.
# most recent common ancestor (MRCA)
nodes_assemblies <- tree_taxonomy %>%
    filter(label %in% c("179_assembly","337_assembly","342_assembly", "543_assembly"))

assemblies_mrca <- MRCA(tree,nodes_assemblies$node )

claude_mrca <- groupClade(as_tibble(tree), assemblies_mrca) %>%
    filter(group==1) %>% 
    left_join(taxonomy_w,
              by=c("label"="gtdb_id")) %>% 
    as.treedata()

write.tree(as.phylo(claude_mrca), file = "claude_mrca.tree", append = FALSE,
           digits = 10, tree.names = FALSE)
# the node 82141 is the family
claude_f__DSM_18226 <- groupClade(as_tibble(tree), 82141) %>%
    filter(group==1) %>% 
    left_join(taxonomy_w,
              by=c("label"="gtdb_id")) %>% 
    as.treedata()

tree_long <- as_tibble(tree)
as.treedata(tree_long)

# visualisation

tree_plot_mrca <- ggtree(data = claude_mrca,
               aes(color=g)) + 
     geom_tiplab(size=2, aes(label= s, color=g)) + 
     theme_tree2(legend.position = c(.1, .88))
    
ggsave(plot=tree_plot_mrca,
       "../tree_plot_mrca.pdf",
       device="pdf",
       height = 200,
       width=150,
       units="cm", limitsize = F)
