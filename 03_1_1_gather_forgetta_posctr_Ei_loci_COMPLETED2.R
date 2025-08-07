#' Obtain loci that contain a gene that is a known positive control for at least 1 of the 12 traits
#' 
#' 1. obtain the positive control gene list and all Ei predictions
#' 2. For each gene with Ei prediction (may be the same gene for different traits with different scores), determine if it is in the positive control gene list
#'    -> not matching by trait, so if a gene is a target for any of the 12 traits, it is labelled as yes
#' 3. Obtain the unique loci for all these genes for each trait
#'  -> Different trait can have similar loci (shifted a bit) and overlap the same gene.

library(data.table)
library(tidyverse)
library(glue)
library(yaml)
library(ggExtra)

args  <- read_yaml("/config/config.yaml")[['proj_resource']]

forgetta_targets <- readRDS("/forgetta_target.Rds")

gene_ei_loci_information  <- read.table(
  "/ei_predictions.txt",sep="\t",quote="",header=T
) %>% dplyr::select(gene_name,trait_name,locus_chromosome,locus_start,locus_end) %>% unique()

gene_ei_loci_information <- gene_ei_loci_information %>% mutate(
  is_forgetta_target = ifelse(
    gene_name %in% forgetta_targets$Gene.Symbol,1,0
  )
)
target_loci  <- gene_ei_loci_information %>% filter(is_forgetta_target == 1) %>%
  dplyr::select(
    trait_name,locus_chromosome,locus_start,locus_end
  ) %>% unique()

saveRDS(
  target_loci,
  "/Forgetta_posCtr_Forgetta_loci_repurposing.Rds"
)
