#!/usr/bin/env Rscript

###############################################################################
# script name: kushneria_rna.R
# developed by: Savvas Paragkamian
# framework: SarrisLab
###############################################################################
# GOAL:
# Aim of this script is to analyse and  visualise the micro and long RNAs 
# of plant primming
###############################################################################
###############################################################################

library(tidyverse)
library(readxl)

mirna_target <- read_delim("../kushneria_rna/miRNA-target_230b9a9a19001.tsv", delim="\t")

mirna <- read_delim("../kushneria_rna/1751277461560_230b9999b3001.tsv", delim="\t")

unique(mirna$`Gene ID`)
unique(mirna_target$`Gene ID`)


not_shared_mi_t <- setdiff(mirna$`Gene ID`,mirna_target$`Gene ID`)
not_shared_t_mi <- setdiff(mirna_target$`Gene ID`,mirna$`Gene ID`)

## summary

mirna_all <- mirna_target |> 
    left_join(mirna, by=c("Gene ID"="Gene ID"))


## mirna target

mirna_distinct <- mirna_all |>
    distinct(`Gene ID`, `miRNA targeted lncRNA or mRNA`)

mirna_target_TAPIR <- mirna_all |>
    filter(!is.na(`TAPIR Score`))

mirna_target_Targetfinder <- mirna_all |>
    filter(!is.na(`Targetfinder Score `))

mirna_all_summaries <- mirna_all |>
    group_by(`Gene ID`) |>
    summarise(n_targets=n()) |>
    arrange(desc(n_targets))

mirna_all_summaries_kegg <- mirna_all |>
    group_by(`Gene ID`,`KEGG Pathway Term`) |>
    summarise(n_targets=n(),.groups="keep") |>
    arrange(desc(n_targets))

mirna_all_summaries_kegg_dist <- mirna_all_summaries_kegg |>
    group_by(`KEGG Pathway Term`) |>
    summarise(n_gene=n(), n_targets=sum(n_targets)) |> 
    arrange(desc(n_targets))

kegg_summary <- mirna_all |>
    group_by(`KEGG Pathway Term`) |>
    summarise(n_mirna=n()) |>
    arrange(desc(n_mirna))




