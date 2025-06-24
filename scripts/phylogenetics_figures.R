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
library(treeio)
library(ggtreeExtra)

# katerina

taxonomy <- read_tsv("../katerina/gtdbtk.bac120.decorated.tree-taxonomy", col_names = c("genome_id", "taxonomy"))
taxonomy_table <- read_tsv("../katerina/gtdbtk.bac120.decorated.tree-table")


family_genomes <- taxonomy %>%
  filter(str_detect(taxonomy, "f__Sphingobacteriaceae")) |>
  pull(genome_id)

genus_genomes <- taxonomy %>%
  filter(str_detect(taxonomy, "g__Olivibacter|g__Pseudosphingobacterium")) |>
  pull(genome_id)

## read the de novo tree
tree <- read.tree("../katerina/gtdbtk.bac120.decorated.tree")

matching_tips <- intersect(tree$tip.label, family_genomes)
Sphingobacteriaceae_tree <- keep.tip(tree, matching_tips)
matching_tips_g <- intersect(tree$tip.label, genus_genomes)
Pseudosphingobacterium_tree <- keep.tip(tree, matching_tips_g)

Sphingobacteriaceae_taxonomy <- taxonomy |>
    filter(genome_id %in% family_genomes) |>
    mutate(genus=str_extract(taxonomy, "g__[^;]+")) |>
    mutate(species=str_extract(taxonomy, "s__[^;]+")) |>
    mutate(species= if_else(is.na(species),genome_id,species))

tree_df <- data.frame(label = Sphingobacteriaceae_tree$tip.label)

# Join with taxonomy
tip_data <- left_join(tree_df, Sphingobacteriaceae_taxonomy, by = c("label" = "genome_id"))

tree_p <- ggtree(Sphingobacteriaceae_tree) %<+% tip_data + 
  geom_tippoint(aes(color = genus), size = 2) +
  geom_tiplab(aes(label = species), size = 1.5) +
  theme(legend.position = "inside",legend.position.inside = c(0.2,0.7))

ggsave(plot=tree_p,
       "../Sphingobacteriaceae_tree_p.pdf",
       device="pdf",
       height = 50,
       width=30,
       units="cm", limitsize = F)
##### genus
matching_tips_g <- intersect(tree$tip.label, genus_genomes)
Pseudosphingobacterium_tree <- keep.tip(tree, matching_tips_g)

genus_taxonomy <- taxonomy |>
    filter(genome_id %in% genus_genomes) |>
    mutate(genus=str_extract(taxonomy, "g__[^;]+")) |>
    mutate(species=str_extract(taxonomy, "s__[^;]+")) |>
    mutate(species= if_else(is.na(species),genome_id,species))

tree_gene_df <- data.frame(label = Pseudosphingobacterium_tree$tip.label)

# Join with taxonomy
tip_data <- left_join(tree_df, genus_taxonomy, by = c("label" = "genome_id"))

xlim_tree <- max(node.depth.edgelength(Pseudosphingobacterium_tree)) + 0.08
tree_p <- ggtree(Pseudosphingobacterium_tree,layout = "rectangular") %<+% tip_data + 
  geom_tippoint(aes(color = genus),size=6) +
  geom_tiplab(
    aes(label = species),
    size = 4,
    align = TRUE,     # extend labels radially outward
    linesize = 0.2,
    offset = 0.01      # push labels slightly outward from tips
  ) +
  coord_cartesian(clip = "off") +  # allows labels outside plot boundary
  xlim(0, xlim_tree) +
  theme(legend.position = "bottom",
        plot.margin = margin(1,1,1,1, "cm")

  ) 

ggsave(plot=tree_p,
       "../Pseudosphingobacterium_tree_p.png",
       device="png",
       height = 15,
       width = 20,
       units="cm", limitsize = F)

# kushneria
#
taxonomy <- read_tsv("../kushneria_denovo/gtdbtk.bac120.decorated.tree-taxonomy", col_names = c("genome_id", "taxonomy"))
taxonomy_table <- read_tsv("../kushneria_denovo/gtdbtk.bac120.decorated.tree-table")

family_genomes <- taxonomy %>%
  filter(str_detect(taxonomy, "f__Halomonadaceae")) |>
  pull(genome_id)

## read the de novo tree
tree <- read.tree("../kushneria_denovo/gtdbtk.bac120.decorated.tree")

matching_tips <- intersect(tree$tip.label, family_genomes)
Halomonadaceae_tree <- keep.tip(tree, matching_tips)

Halomonadaceae_taxonomy <- taxonomy |>
    filter(genome_id %in% family_genomes) |>
    mutate(genus=str_extract(taxonomy, "g__[^;]+")) |>
    mutate(species=str_extract(taxonomy, "s__[^;]+")) |>
    mutate(species= if_else(is.na(species),genome_id,species))

tree_df <- data.frame(label = Halomonadaceae_tree$tip.label)

# Join with taxonomy
tip_data <- left_join(tree_df, Halomonadaceae_taxonomy, by = c("label" = "genome_id"))

Halomonadaceae_tree_p <- ggtree(Halomonadaceae_tree) %<+% tip_data + 
  geom_tippoint(aes(color = genus), size = 2) +
  geom_tiplab(aes(label = species), size = 1.5) +
  theme(legend.position = "inside",legend.position.inside = c(0.2,0.7))

ggsave(plot=Halomonadaceae_tree_p,
       "../Halomonadaceae_tree_p.pdf",
       device="pdf",
       height = 50,
       width=30,
       units="cm", limitsize = F)

## Kushneria only


ggsave(plot=Halomonadaceae_tree_p,
       "../Halomonadaceae_tree_p.png",
       device="png",
       dpi=300,
       height = 60,
       width=40,
       units="cm", limitsize = F)

# bacillus 
# read the taxonomy

# gawk -F"\t" '{taxon=$1; split($2,taxonomy, "; ") ; print "gtdb_id" "\t" "classification_full" "\t" "classification" "\t" "taxonname" ; for (i in taxonomy){split(taxonomy[i], classification,"__"); print taxon "\t" taxonomy[i] "\t" classification[1] "\t" classification[2]}}' gtdbtk.bac120.decorated.tree-taxonomy > gtdbtk.taxonomy_long.tsv

taxonomy <- read_delim( "../bac3new_de_novo/infer/gtdbtk.taxonomy_long.tsv",
                       delim="\t",
                       col_names=T)
# transform taxonomy to wide format
taxonomy_w <- taxonomy %>% group_by(gtdb_id, classification) %>%
    summarise(taxa=paste(taxonname, collapse = "; "),
              .groups = "drop") %>%  
    pivot_wider(names_from=classification, values_from=taxa)

# load tree
tree <- ape::read.tree(file = "../bac3new_de_novo/gtdbtk.bac120.decorated.tree")
tree_table <- read_delim( "../bac3new_de_novo/gtdbtk.bac120.decorated.tree-table", delim="\t")
tree_taxonomy <- as_tibble(tree) %>%
    left_join(taxonomy_w, by=c("label"="gtdb_id")) %>%
    treeio::as.treedata()

# filter the tree

tree_bacilli <- as_tibble(tree) %>% filter(grepl("Bacilli",label))

# it doesn't have a node for Bacilli.
# most recent common ancestor (MRCA)
nodes_assemblies <- as_tibble(tree) %>%
    filter(label %in% c("179_assembly", "543_assembly")) #"337_assembly"

all_assemblies <- as_tibble(tree) %>%
    filter(label %in% c("179_assembly", "543_assembly","337_assembly")) #"337_assembly"

all_assemblies_m <- MRCA(tree,all_assemblies$node ) 
assemblies_mrca <- MRCA(tree,nodes_assemblies$node )

tree_assemblies <- extract.clade(tree, assemblies_mrca) %>%
    as_tibble() %>%
    left_join(taxonomy_w, by=c("label"="gtdb_id")) %>%
    as.treedata()

all_assemblies <- extract.clade(tree, all_assemblies_m) %>%
    as_tibble() %>%
    left_join(taxonomy_w, by=c("label"="gtdb_id")) %>%
  mutate(species=ifelse(s=="NA",label,s)) %>%
  as.treedata()

#claude_mrca <- groupClade(tree_taxonomy, assemblies_mrca) %>%
#    filter(group==1) 

#write.tree(as.phylo(claude_mrca), file = "claude_mrca.tree", append = FALSE,
#           digits = 10, tree.names = FALSE)

# the node 82141 is the family
claude_f__DSM_18226 <- groupClade(as_tibble(tree), 82141) %>%
    filter(group==1) %>% 
    left_join(taxonomy_w,
              by=c("label"="gtdb_id")) %>% 
    as.treedata()

# visualisation

tree_plot_mrca <- ggtree(tree_assemblies) + 
     geom_tiplab(size=2, aes(as_ylab=TRUE, color=g)) + 
#     geom_text(aes(label=branch.length, color=g)) + 
     theme_tree2(legend.position = c(.1, .88))
    
ggsave(plot=tree_plot_mrca,
       "../tree_plot_mrca_gtdb.pdf",
       device="pdf",
       height = 80,
       width=40,
       units="cm", limitsize = F)

tree_plot_mrca2 <- ggtree(tree_assemblies) + 
     geom_tiplab(size=2, aes(label=s, color=g)) + 
#     geom_text(aes(label=branch.length, color=g)) + 
     theme_tree2(legend.position = c(.1, .88))

ggsave(plot=tree_plot_mrca2,
       "../tree_plot_mrca_s.pdf",
       device="pdf",
       height = 80,
       width=40,
       units="cm", limitsize = F)


tree_plot_mrca3 <- ggtree(all_assemblies) + 
     geom_tiplab(size=2, aes(label=species, color=g)) + 
#     geom_text(aes(label=branch.length)) + 
     theme_tree2(legend.position = c(.1, .88))

#tree_plot_mrca3 <- ggtree(all_assemblies) + 
#  geom_tiplab(size=2, aes(label=species)) + 
#  #     geom_text(aes(label=branch.length)) + 
#  theme_tree2(legend.position = c(.1, .88))

ggsave(plot=tree_plot_mrca3,
       "../tree_plot_mrca_3.pdf",
       device="pdf",
       height = 100,
       width=45,
       units="cm", limitsize = F)

