#' for the relevant trait
#' 1. Pick the genes with largest Ei scores per trait, locus (locus chr, start, end)
#' 2. filter to be in proteome gene list
#' 

library(tidyverse)
library(glue)
library(data.table)
library(yaml)

args  <- read_yaml("/config/config.yaml")[['proj_resource']]


outdir  <- file.path('/ei_genes')




ei_trait_name  <- args$ei_12_traits_map[args$ei_12_traits]

ei_data  <- read.table(
  "/ei_predictions.txt",sep="\t",header=T,quote=""
)
assertthat::assert_that(
  all(
    ei_trait_name %in% unique(ei_data$trait_name)
  )
)

gene_annotation_file  <- readRDS("/relevant_ensembl_canonical_protein_coding.Rds")


for (each_trait in ei_trait_name){
  print(each_trait)
  ei_trait_predictions  <- ei_data %>% filter(trait_name == each_trait)
  best_gene  <- ei_trait_predictions %>% group_by(locus_chromosome,locus_start,locus_end) %>% slice_max(Ei_prediction) %>% ungroup()
  best_gene  <- best_gene %>% filter(
    gene_name %in% gene_annotation_file$gene_name
  )
  saveRDS(
    best_gene,
    file.path(
      outdir,glue("{each_trait}_best_ei_gene.Rds")
    )
  )
}